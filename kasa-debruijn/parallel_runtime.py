"""Bounded Windows-spawn processes. This module contains no construction math.

Only small task descriptions cross queues. Large buffers use file mappings.
On failure/interrupt, owned workers are terminated and joined before cleanup.
"""
import multiprocessing as mp
import os
import pickle
import queue
import signal
import traceback
import time

_progress_sink=None
_progress_last=0.0

def report_progress(done,total):
    global _progress_last
    if _progress_sink is None:return
    now=time.monotonic()
    if now-_progress_last>=0.5:
        _progress_last=now; _progress_sink(done,total)


def process_memory(pid):
    """Windows working set and process-lifetime peak; own processes only."""
    if os.name != 'nt': return None
    import ctypes
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_=[('cb',wintypes.DWORD),('faults',wintypes.DWORD)]+[(name,ctypes.c_size_t) for name in ('peak','working','quota_peak_paged','quota_paged','quota_peak_nonpaged','quota_nonpaged','pagefile','peak_pagefile')]
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    psapi=ctypes.WinDLL('psapi',use_last_error=True)
    kernel.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD]; kernel.OpenProcess.restype=wintypes.HANDLE
    kernel.CloseHandle.argtypes=[wintypes.HANDLE]
    psapi.GetProcessMemoryInfo.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
    handle=kernel.OpenProcess(0x410,False,pid)
    if not handle: return None
    try:
        counters=Counters(); counters.cb=ctypes.sizeof(counters)
        if not psapi.GetProcessMemoryInfo(handle,ctypes.byref(counters),counters.cb): return None
        return counters.working,counters.peak
    finally: kernel.CloseHandle(handle)


def worker_limit(requested, work, tasks, memory_limit_mib=8192, per_worker_bytes=64*1024**2, parent_bytes=128*1024**2):
    if type(requested) is not int or requested < 0 or requested > 64:
        raise ValueError('--workers must be 0 (automatic) or 1..64')
    if type(memory_limit_mib) is not int or memory_limit_mib<64: raise ValueError('memory limit must be at least 64 MiB')
    count=(4 if work >= 8_000_000 else 1) if requested == 0 else requested
    budget=memory_limit_mib*1024**2
    memory_workers=max(1,(budget-parent_bytes)//per_worker_bytes)
    return max(1,min(count,os.cpu_count() or 1,tasks,memory_workers))


def _worker(inputs,outputs):
    signal.signal(signal.SIGINT,signal.SIG_IGN)
    while True:
        task=inputs.get()
        if task is None: return
        batch,index,fn,args,active=task
        global _progress_sink,_progress_last
        _progress_last=0.0
        _progress_sink=(lambda done,total:outputs.put(('progress',batch,index,done,total)))if active else None
        try:
            result=fn(*args)
            pickle.dumps(result)
        except BaseException:
            outputs.put((batch,index,False,traceback.format_exc()))
            return
        outputs.put((batch,index,True,result))


class WorkerPool:
    def __init__(self,workers,memory_limit_mib=8192,temp_dir=None):
        self.workers=workers
        self.temp_dir=temp_dir
        self.processes=[]; self.batch=0; self.closed=False
        self.inputs=None; self.outputs=None
        self.memory_limit=memory_limit_mib*1024**2
        self.peak_worker_sum=0; self.peak_sampled_sum=0

    def sample_memory(self):
        parent=process_memory(os.getpid()); children=[process_memory(p.pid) for p in self.processes]
        if parent is None:
            if os.name=='nt': raise RuntimeError('cannot measure this process working set for the memory budget')
            return
        if any(x is None for x in children):
            if any(x is None and p.exitcode is None for p,x in zip(self.processes,children)):
                raise RuntimeError('cannot measure an active worker working set for the memory budget')
            return
        self.peak_worker_sum=max(self.peak_worker_sum,sum(x[1] for x in children))
        self.peak_sampled_sum=max(self.peak_sampled_sum,parent[0]+sum(x[0] for x in children))
        if parent[1]+self.peak_worker_sum>self.memory_limit:
            raise RuntimeError('memory budget exceeded: measured %d bytes, allowed %d bytes; raise --memory-limit-mib'%(parent[1]+self.peak_worker_sum,self.memory_limit))

    def memory_stats(self):
        parent=process_memory(os.getpid())
        total=None if parent is None else parent[1]+self.peak_worker_sum
        if total is not None and total>self.memory_limit: raise RuntimeError('memory budget exceeded: measured %d bytes, allowed %d bytes; raise --memory-limit-mib'%(total,self.memory_limit))
        return dict(limit_mib=self.memory_limit//1024**2,parent_process_lifetime_peak_bytes=None if parent is None else parent[1],worker_peak_working_set_sum_bytes=self.peak_worker_sum,conservative_peak_total_bytes=total,sampled_simultaneous_peak_bytes=self.peak_sampled_sum,method='Windows process working sets; shared mapped pages may be counted more than once; sum of per-process peaks is conservative')

    def _start(self):
        context=mp.get_context('spawn')
        self.inputs=context.Queue(); self.outputs=context.Queue()
        try:
            for _ in range(self.workers):
                process=context.Process(target=_worker,args=(self.inputs,self.outputs))
                process.start(); self.processes.append(process)
        except BaseException:
            self.close(abort=True)
            raise

    def map(self,fn,tasks,progress=None,task_units=None):
        global _progress_sink,_progress_last
        if self.closed: raise RuntimeError('worker pool is closed')
        tasks=list(tasks)
        for args in tasks:
            if len(pickle.dumps((fn,args)))>4*1024**2: raise ValueError('worker task exceeds bounded metadata size')
        units=list(task_units)if task_units is not None else [1]*len(tasks)
        if len(units)!=len(tasks) or any(x<=0 for x in units):raise ValueError('invalid progress units')
        amounts=[0]*len(tasks); completed=set(); total=sum(units)
        def notify(index,done):
            if not 0<=done<=units[index]:raise RuntimeError('invalid worker progress count')
            amounts[index]=max(amounts[index],done)
            if progress:progress(len(completed),len(tasks),sum(amounts),total,index)
        if progress:progress(0,len(tasks),0,total,None)
        if self.workers == 1:
            out=[]
            try:
                for index,args in enumerate(tasks):
                    _progress_last=0.0
                    _progress_sink=(lambda done,n,i=index:notify(i,done))if progress else None
                    out.append(fn(*args)); completed.add(index);notify(index,units[index])
            finally:_progress_sink=None
            self.sample_memory(); return out
        if not self.processes: self._start()
        self.batch+=1; batch=self.batch
        for index,args in enumerate(tasks): self.inputs.put((batch,index,fn,args,progress is not None))
        results=[None]*len(tasks); received=set()
        try:
            while len(received) != len(tasks):
                self.sample_memory()
                try:
                    message=self.outputs.get(timeout=0.1)
                    if message[0]=='progress':
                        _,b,index,done,n=message
                        if b!=batch or not 0<=index<len(tasks) or n!=units[index]:raise RuntimeError('invalid worker progress message')
                        notify(index,done);continue
                    b,index,ok,value=message
                except queue.Empty:
                    failed=[p.exitcode for p in self.processes if p.exitcode is not None]
                    if failed: raise RuntimeError('worker exited before completing its batch: '+str(failed))
                    if progress:progress(len(completed),len(tasks),sum(amounts),total,None)
                    continue
                if b != batch or index in received or not 0 <= index < len(tasks):
                    raise RuntimeError('invalid or duplicate worker result')
                if not ok: raise RuntimeError('worker failed:\n'+value)
                results[index]=value; received.add(index); completed.add(index);notify(index,units[index])
            self.sample_memory()
        except BaseException:
            self.close(abort=True)
            raise
        return results

    def close(self,abort=False):
        if self.closed: return
        self.closed=True
        if self.processes:
            if abort:
                for p in self.processes:
                    if p.is_alive(): p.terminate()
            else:
                for _ in self.processes: self.inputs.put(None)
            for p in self.processes:
                p.join(timeout=5)
                if p.is_alive(): p.kill(); p.join()
                p.close()
        for q in (self.inputs,self.outputs):
            if q is not None:
                q.cancel_join_thread(); q.close()

    def __enter__(self): return self
    def __exit__(self,kind,value,tb): self.close(abort=kind is not None)
