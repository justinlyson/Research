"""Deterministic revised construction and atomic, process-parallel family export."""
import argparse
import hashlib
import json
import mmap
import os
from pathlib import Path
import platform
import tempfile
import time
from construction import construct,require
from parallel_runtime import WorkerPool,worker_limit,report_progress
from progress import Progress


def export_phase(reference,directory,q,N,phase):
    path=Path(directory)/('cycle_%03d.txt' % phase)
    table=[str(0 if x == 0 else 1+(x-1+phase) % (q-1)) for x in range(q)]
    h=hashlib.sha256()
    with open(reference,'rb') as source,mmap.mmap(source.fileno(),0,access=mmap.ACCESS_READ) as a,path.open('xb') as f:
        for i in range(0,N,4096):
            data=(' '.join(table[x] for x in a[i:i+4096])+'\n').encode('ascii')
            f.write(data); h.update(data)
            report_progress(min(i+4096,N-1),N)
    return dict(name=path.name,text_sha256=h.hexdigest(),bytes=path.stat().st_size)


def generate(q,k,directory,max_states=100_000_000,workers=0,memory_limit_mib=8192,temp_dir=None,quiet=False):
    total_start=time.perf_counter(); display=Progress(quiet)
    directory=Path(directory).absolute()
    require(not directory.exists(),'output directory already exists; choose a new path')
    require(type(q) is int and 3 <= q <= 255 and type(k) is int and k >= 2,'revised construction requires 3<=q<=255 and k>=2')
    require(type(max_states) is int and max_states>0,'--max-states must be a positive integer')
    require(q**k<=0xffffffff,'uint32 implementation bound exceeded: q^k must be <= 4294967295; this is not a removable resource budget')
    require(q**k<=max_states,'state budget exceeded: requires %d states, allowed %d; raise --max-states'%(q**k,max_states))
    if temp_dir is not None:
        temp_dir=Path(temp_dir)
        require(temp_dir.is_dir(),'--temp-dir must name an existing directory: '+str(temp_dir))
    N=q**k; count=worker_limit(workers,N*(q-1),q-1,memory_limit_mib,64*1024**2+N,128*1024**2+20*N)
    directory.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.kasa-build-',dir=directory.parent) as tmp:
        staging=Path(tmp)/'complete'; staging.mkdir()
        with WorkerPool(count,memory_limit_mib,temp_dir) as pool:
            start=time.perf_counter(); P,construction=construct(q,k,max_states,pool,display)
            construction_seconds=time.perf_counter()-start
            start=time.perf_counter(); prep=display.stage('Reference cycle preparation',N,'positions','writing'); a=bytearray(N); v=0; head_mod=q**(k-1)
            for i in range(N):
                x=v//head_mod
                a[i]=x if q == 3 else 0 if x == q-1 else x+1
                v=P[v]
                if i & 0x3ffff==0:prep.update(i)
            require(v == 0,'reference cycle did not close')
            pos=(a+a[:k-1]).find(bytes(k))
            require(pos >= 0,'zero-word root missing')
            a=a[pos:]+a[:pos]
            reference=Path(tmp)/'reference.u8'; reference.write_bytes(a)
            reference_preparation_seconds=time.perf_counter()-start
            del P,a
            prep.finish()
            start=time.perf_counter(); export=display.stage('Export phases 0..%d'%(q-2),N*(q-1),'symbols','writing')
            files=pool.map(export_phase,[(str(reference),str(staging),q,N,phase) for phase in range(q-1)],progress=export.pool_update if display.enabled else None,task_units=[N]*(q-1))
            writing_seconds=time.perf_counter()-start+reference_preparation_seconds
            memory=pool.memory_stats()
        export.finish('%d/%d phase files written'%(q-1,q-1))
        require([f['name'] for f in files] == ['cycle_%03d.txt' % t for t in range(q-1)],'worker export coverage differs')
        result=dict(revision='runtime-2026.10',q=q,k=k,format='q-1 numeric-token cyclic words; exactly q^k symbols; literal 0 fixed',
                    temp_directory=None if temp_dir is None else str(temp_dir),
                    construction=construction,construction_seconds=construction_seconds,writing_seconds=writing_seconds,
                    reference_preparation_seconds=reference_preparation_seconds,workers_requested=workers,workers_used=count,
                    memory=memory,
                    construction_parallel_threshold_states=2_000_000,
                    machine=dict(platform=platform.platform(),python=platform.python_version(),processor=platform.processor()),files=files)
        (staging/'generation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        # Windows rename fails if another creator wins; it never overwrites it.
        if os.name == 'nt':
            os.rename(staging,directory)
        else:
            directory.mkdir(exist_ok=False); committed=[]
            try:
                for f in staging.iterdir():
                    target=directory/f.name
                    with target.open('xb') as out,f.open('rb') as source:
                        committed.append(target)
                        while data:=source.read(1024*1024): out.write(data)
            except BaseException:
                for f in committed: f.unlink(missing_ok=True)
                try: directory.rmdir()
                except OSError: pass
                raise
    result['total_generation_seconds']=time.perf_counter()-total_start
    result['overhead_seconds']=max(0,result['total_generation_seconds']-construction_seconds-writing_seconds)
    result['stage_timings']=display.timings
    (directory/'generation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    display.log('Generation complete: %d/%d cycles; total %.3fs (construction %.3fs; writing/preparation %.3fs)'%(q-1,q-1,result['total_generation_seconds'],construction_seconds,writing_seconds))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--q',type=int,required=True); p.add_argument('--k',type=int,required=True); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--max-states',type=int,default=100_000_000)
    p.add_argument('--workers',type=int,default=0,help='0=automatic (small jobs 1, large jobs up to 4); 1=serial; maximum 64')
    p.add_argument('--memory-limit-mib',type=int,default=8192)
    p.add_argument('--temp-dir',type=Path,help='existing directory for raw-lift temporary mappings; export staging stays beside --out')
    p.add_argument('--quiet',action='store_true',help='suppress live stderr progress')
    p.add_argument('--json',action='store_true',help='print the full generation receipt')
    a=p.parse_args()
    try: result=generate(a.q,a.k,a.out,a.max_states,a.workers,a.memory_limit_mib,a.temp_dir,a.quiet)
    except (ValueError,OSError,RuntimeError) as e: p.exit(1,'FAIL: '+str(e)+'\n')
    except KeyboardInterrupt: p.exit(130,'Interrupted. Any completed atomic export is retained; verify it before use.\n')
    if a.json: print(json.dumps(result,indent=2))
    else:
        print('Generated %d cycles for q=%d,k=%d using %d worker(s).'%(a.q-1,a.q,a.k,result['workers_used']))
        print('Total generation %.3f s.'%result['total_generation_seconds'])
        print('Construction %.3f s; export preparation/writing %.3f s.'%(result['construction_seconds'],result['writing_seconds']))
        print('Output: '+str(a.out))


if __name__ == '__main__': main()
