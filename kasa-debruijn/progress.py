"""Throttled parent-owned stderr progress and measured stage timings."""
import sys
import time


class Progress:
    def __init__(self,quiet=False):
        self.enabled=not quiet; self.timings=[]

    def stage(self,label,total,unit='items',group='construction'):
        return Stage(self,label,total,unit,group)

    def log(self,text):
        if self.enabled: print(text,file=sys.stderr,flush=True)


class Stage:
    def __init__(self,owner,label,total,unit,group):
        self.owner=owner; self.label=label; self.total=max(1,total)
        self.unit=unit; self.group=group; self.start=time.perf_counter(); self.last=0
        self.update(0,force=True)

    def update(self,done,detail='',force=False):
        if not self.owner.enabled:return
        now=time.perf_counter()
        if not force and now-self.last<1:return
        self.last=now; elapsed=now-self.start; done=min(self.total,max(0,done))
        percent=100*done/self.total
        percent_text=('%.1f%%'%percent)if percent<99.9 else '>99.9% (work remaining)'
        if done==self.total:percent_text='finalizing'
        rate=done/elapsed if elapsed else 0
        eta=('; ETA ~%.1fs'%((self.total-done)/rate)) if rate and done<self.total else ''
        speed=('; %.0f %s/s'%(rate,self.unit))if rate else ''
        print('%s: %s (%s/%s %s); elapsed %.1fs%s%s%s'%(self.label,percent_text,format(done,','),format(self.total,','),self.unit,elapsed,speed,eta,('; '+detail)if detail else ''),file=sys.stderr,flush=True)

    def pool_update(self,done_tasks,tasks,done,total,index=None):
        self.update(done,'%d/%d cycles/chunks complete%s'%(done_tasks,tasks,('; phase/chunk %d'%index)if index is not None else ''))

    def finish(self,detail=''):
        seconds=time.perf_counter()-self.start
        self.owner.timings.append(dict(stage=self.label,group=self.group,seconds=seconds))
        if self.owner.enabled:
            print('%s: 100.0%% (%s/%s %s); elapsed %.3fs%s'%(self.label,format(self.total,','),format(self.total,','),self.unit,seconds,('; '+detail)if detail else ''),file=sys.stderr,flush=True)
