"""Independent multicore verifier. Imports no construction or generator code.

Every file is checked; the parent explicitly checks the global arc union.
q/k are caller arguments, never taken from generator metadata.
"""
import argparse
from collections import Counter
import hashlib
import json
import mmap
from pathlib import Path
import tempfile
import time
from parallel_runtime import WorkerPool,worker_limit,report_progress
from progress import Progress


def require(ok,message):
    if not ok: raise ValueError(message)


def read_cycle(path,q,N):
    symbols=bytearray()
    with Path(path).open('rb') as f:
        for line in f:
            if len(line)>65536: raise ValueError('line too long')
            tokens=line.split()
            if len(symbols)+len(tokens)>N: raise ValueError('too many cyclic positions')
            for token in tokens:
                if not token.isdigit(): raise ValueError('noninteger alphabet token')
                x=int(token)
                if x>=q: raise ValueError('symbol outside alphabet')
                symbols.append(x)
    if len(symbols)!=N: raise ValueError('incorrect number of cyclic positions')
    return symbols


def _check_cycle(path,q,k,N,phase,reference,arcfile,stride,symmetry):
    a=read_cycle(path,q,N)
    digest=hashlib.sha256(a).hexdigest()
    if symmetry:
        table=bytes(0 if x==0 else 1+(x-1+phase) % (q-1) for x in range(256))
        with open(reference,'rb') as f,mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as ref:
            for start in range(0,N,65536):
                if a[start:start+65536] != ref[start:start+65536].translate(table):
                    raise ValueError('prescribed symmetry differs at phase '+str(phase))
    seen=bytearray(N); arcs=bytearray(N*stride); v=0
    for x in a[:k]: v=v*q+x
    root=v; mod=q**(k-1); a.extend(a[:k]); next_letters=memoryview(a)[k:]
    try:
        for start in range(0,N,262144):
            for y in next_letters[start:start+262144]:
                if seen[v]: raise ValueError('repeated length-k cyclic word in '+Path(path).name)
                seen[v]=1
                index=v*stride+(y>>3)
                arcs[index] |= 1 << (y&7)
                v=(v % mod)*q+y
            report_progress(min(start+262144,N-1),N)
    finally: next_letters.release()
    if v!=root: raise ValueError('cyclic closure failed')
    with open(arcfile,'xb') as f: f.write(arcs)
    return dict(file=Path(path).name,positions=N,all_words_once=True,wraparound_checked=True,symbols_sha256=digest,arcfile=str(arcfile),phase=phase)


def residual_report(packed,q,k,N,stride,stage=None):
    mod=q**(k-1); full=(1<<q)-1
    residual=bytearray(N); incoming=bytearray(N)
    format_code={1:'B',2:'H',4:'I',8:'Q'}.get(stride)
    view=memoryview(packed).cast(format_code) if format_code else None
    try:
        for v in range(N):
            if stage is not None and v & 0x3ffff==0:stage.update(v,'checking missing arcs and indegrees')
            mask=view[v] if view is not None else int.from_bytes(packed[v*stride:(v+1)*stride],'little')
            if mask & ~full: raise ValueError('global arc bitset uses an invalid alphabet bit')
            missing=full ^ mask
            if not missing or missing & (missing-1): raise ValueError('residual outdegree differs from 1')
            y=missing.bit_length()-1; residual[v]=y
            target=(v % mod)*q+y
            if incoming[target]: raise ValueError('residual indegree differs from 1')
            incoming[target]=1
    finally:
        if view is not None: view.release()
    constants={x*((N-1)//(q-1)) for x in range(q)}
    for v in constants:
        if (v % mod)*q+residual[v]!=v: raise ValueError('forced constant loop absent')
    visited=bytearray(N); histogram=Counter(); loops=[]; representatives=[]
    def decode(v):
        w=[0]*k
        for i in range(k-1,-1,-1): v,w[i]=divmod(v,q)
        return w
    visited_total=0
    for root in range(N):
        if visited[root]: continue
        v=root; length=0
        while not visited[v]:
            visited[v]=1; length+=1; v=(v % mod)*q+residual[v]
            if stage is not None and length & 0x3ffff==0:stage.update(N+visited_total+length,'traversing residual cycles')
        visited_total+=length
        if stage is not None:stage.update(N+visited_total,'traversing residual cycles')
        if v!=root: raise ValueError('residual component does not close at root')
        if root in constants:
            if length!=1: raise ValueError('constant component is not a loop')
            loops.append(root)
        else:
            histogram[length]+=1
            representatives.append(dict(length=length,representative_word=decode(root)))
    summary=dict(indegree=1,outdegree=1,constant_loops=len(loops),constant_loop_words=[decode(v) for v in sorted(loops)],nonconstant_cycle_count=sum(histogram.values()),nonconstant_cycle_lengths=dict(sorted(histogram.items())),nonconstant_vertices=sum(l*c for l,c in histogram.items()))
    return summary,representatives


def resource_plan(q,k,max_states=100_000_000,max_arcs=10_000_000_000,max_temp_mib=16384):
    require(3<=q<=255 and k>=2,'verifier domain: 3<=q<=255,k>=2')
    for name,value in (('--max-states',max_states),('--max-arcs',max_arcs),('--max-temp-mib',max_temp_mib)):
        require(type(value) is int and value>0,name+' must be a positive integer')
    N=q**k
    require(N<=0xffffffff,'uint32 implementation bound exceeded: q^k must be <= 4294967295; this is not a removable resource budget')
    require(N<=max_states,'state budget exceeded: requires %d states, allowed %d; raise --max-states'%(N,max_states))
    require(N*q<=max_arcs,'arc budget exceeded: requires %d possible arcs, allowed %d; raise --max-arcs'%(N*q,max_arcs))
    stride=1 if q<=8 else 2 if q<=16 else 4 if q<=32 else 8 if q<=64 else (q+7)//8
    needed=N*stride*(q-1); allowed=max_temp_mib*1024**2; minimum=(needed+1024**2-1)//1024**2
    require(needed<=allowed,'temporary arc bitsets require %d bytes (%d MiB rounded up), allowed %d bytes (%d MiB); set --max-temp-mib to at least %d and ensure free disk space'%(needed,minimum,allowed,max_temp_mib,minimum))
    return dict(states=N,possible_arcs=N*q,stride=stride,temporary_arc_bytes=needed,temporary_reference_bytes=N,minimum_arc_temp_mib=minimum,max_temp_mib=max_temp_mib)


def verify(directory,q,k,symmetry=True,max_states=100_000_000,max_arcs=10_000_000_000,workers=0,memory_limit_mib=8192,max_temp_mib=16384,temp_dir=None,quiet=False):
    started=time.perf_counter(); display=Progress(quiet)
    plan=resource_plan(q,k,max_states,max_arcs,max_temp_mib)
    N=plan['states']; stride=plan['stride']
    if temp_dir is not None:
        temp_dir=Path(temp_dir)
        require(temp_dir.is_dir(),'--temp-dir must name an existing directory: '+str(temp_dir))
    paths=sorted(Path(directory).glob('cycle_*.txt'))
    require(len(paths)==q-1,'export must contain exactly q-1 cycles')
    count=worker_limit(workers,N*(q-1),q-1,memory_limit_mib,32*1024**2+N*(stride+3),64*1024**2+N*(stride*3+6))
    with tempfile.TemporaryDirectory(prefix='kasa-independent-verify-',dir=temp_dir) as tmp:
        root=Path(tmp); reference=root/'reference.u8'
        # This reference comes from the actual text export, never generator state.
        if symmetry: reference.write_bytes(read_cycle(paths[0],q,N))
        checks=display.stage('Per-cycle checks (phases 0..%d)'%(q-2),N*(q-1),'positions','verification')
        with WorkerPool(count,memory_limit_mib) as pool:
            tasks=[(str(path),q,k,N,phase,str(reference),str(root/('arcs_%03d.bin'%phase)),stride,symmetry) for phase,path in enumerate(paths)]
            reports=pool.map(_check_cycle,tasks,progress=checks.pool_update if display.enabled else None,task_units=[N]*(q-1))
            memory=pool.memory_stats()
        checks.finish('%d/%d cycles checked'%(q-1,q-1))
        file_check_seconds=time.perf_counter()-started
        require([r['phase'] for r in reports]==list(range(q-1)),'worker phase coverage differs')
        # Parent-owned exact global union: native integer bit operations are
        # simply a packed set representation, with one bit per (vertex,symbol).
        # This catches collisions BETWEEN workers even with symmetry disabled.
        start=time.perf_counter(); merge=display.stage('GLOBAL collision merge',q-1,'cycle arc sets','verification'); union=0
        for phase,r in enumerate(reports):
            arcfile=root/('arcs_%03d.bin'%phase)
            require(r['arcfile']==str(arcfile) and arcfile.stat().st_size==N*stride,'worker arc file shape differs')
            bits=int.from_bytes(arcfile.read_bytes(),'little')
            require(bits.bit_count()==N,'worker arc coverage differs')
            if union & bits: raise ValueError('repeated length-(k+1) arc between exported cycles (global parent check)')
            union |= bits
            merge.update(phase+1,'global union checked through this set')
            del bits
            arcfile.unlink()
            del r['arcfile']; del r['phase']
        require(union.bit_count()==N*(q-1),'global arc union count differs')
        packed=union.to_bytes(N*stride,'little'); del union
        global_union_seconds=time.perf_counter()-start
        merge.finish()
        start=time.perf_counter(); analysis=display.stage('Residual analysis',2*N,'vertices','verification'); residual,representatives=residual_report(packed,q,k,N,stride,analysis)
        residual_seconds=time.perf_counter()-start
        analysis.finish('%d constant loops; %d other cycles'%(residual['constant_loops'],residual['nonconstant_cycle_count']))
        memory=pool.memory_stats()
    total_seconds=time.perf_counter()-started
    display.log('Verification complete: all %d cycles and GLOBAL arc union passed; total %.3fs'%(q-1,total_seconds))
    return dict(result='PASS',revision='runtime-2026.10',q=q,k=k,cycles_checked=len(paths),positions_per_cycle=N,all_exported_cycles_checked=True,prescribed_symmetry_checked=symmetry,
                resource_plan=plan,temp_directory=None if temp_dir is None else str(temp_dir),
                global_cross_family_union_checked_by_parent=True,workers_requested=workers,workers_used=count,
                memory=memory,
                residual=residual,nonconstant_residual_cycles=representatives,cycles=reports,
                file_check_seconds=file_check_seconds,global_union_seconds=global_union_seconds,residual_seconds=residual_seconds,validation_seconds=total_seconds,stage_timings=display.timings)


def save_report(path,text):
    tmp=None
    try:
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=path.parent,delete=False,prefix='.kasa-report-') as f:
            tmp=Path(f.name); f.write(text+'\n')
        tmp.replace(path)
    finally:
        if tmp is not None: tmp.unlink(missing_ok=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path); p.add_argument('--q',type=int,required=True); p.add_argument('--k',type=int,required=True)
    p.add_argument('--report',type=Path); p.add_argument('--no-symmetry',action='store_true')
    p.add_argument('--max-states',type=int,default=100_000_000); p.add_argument('--max-arcs',type=int,default=10_000_000_000)
    p.add_argument('--workers',type=int,default=0,help='0=automatic; 1=serial fallback; maximum 64')
    p.add_argument('--memory-limit-mib',type=int,default=8192)
    p.add_argument('--max-temp-mib',type=int,default=16384,help='temporary arc-bitset disk budget in MiB; default 16384')
    p.add_argument('--temp-dir',type=Path,help='existing directory for temporary arc bitsets and reference')
    p.add_argument('--quiet',action='store_true',help='suppress live stderr progress')
    p.add_argument('--json',action='store_true',help='print the full result including every residual representative')
    a=p.parse_args()
    try:
        result=verify(a.directory,a.q,a.k,not a.no_symmetry,a.max_states,a.max_arcs,a.workers,a.memory_limit_mib,a.max_temp_mib,a.temp_dir,a.quiet)
        text=json.dumps(result,indent=2)
        if a.report: save_report(a.report,text)
    except (ValueError,OSError,RuntimeError) as e: p.exit(1,'FAIL: '+str(e)+'\n')
    except KeyboardInterrupt: p.exit(130,'Cancelled; verification did not complete.\n')
    if a.json: print(text)
    else:
        print('PASS: every exported cycle and the global arc union checked (%d cycles, %d positions each).'%(result['cycles_checked'],result['positions_per_cycle']))
        print('Workers %d; validation %.3f s.'%(result['workers_used'],result['validation_seconds']))
        print('Residual: %d constant loops; %d other cycles.'%(result['residual']['constant_loops'],result['residual']['nonconstant_cycle_count']))
        if a.report: print('Full report: '+str(a.report))


if __name__ == '__main__': main()
