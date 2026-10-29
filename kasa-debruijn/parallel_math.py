"""Parallel evaluation of the SAME raw Latin lift, in disjoint vertex ranges."""
from array import array
import mmap
from pathlib import Path
import tempfile
from parallel_runtime import report_progress


def _lift_range(pfile,cfile,dfile,ufile,vfile,q,R,D,start,stop):
    handles=[]; maps=[]; views=[]
    try:
        for path,mode in ((pfile,'rb'),(cfile,'rb'),(dfile,'rb'),(ufile,'r+b'),(vfile,'r+b')):
            f=open(path,mode); handles.append(f)
            m=mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ if mode=='rb' else mmap.ACCESS_WRITE)
            maps.append(m); views.append(memoryview(m).cast('I'))
        P,C,previous,U,V=views
        lower_suffix=len(P)//q
        for prefix in range(start,stop):
            d0=previous[prefix]*q; colors=D[prefix % q]
            base=(prefix % lower_suffix)*q*q
            offset=prefix*q
            for last in range(q):
                d=d0+colors[last]; tail=base+last*q
                U[offset+last]=tail+R[P[d] % q][last]
                V[offset+last]=tail+R[C[d] % q][last]
            if (prefix-start) & 4095==0:report_progress(min((prefix-start+1)*q,(stop-start)*q-1),(stop-start)*q)
        return start,stop
    finally:
        # Named aliases also hold exported views until released explicitly.
        for view in views: view.release()
        for m in maps: m.close()
        for f in handles: f.close()


def parallel_lift(P,C,q,k,R,D,pool,stage=None):
    # Prefix projection is only q^k entries; never build a full upper projection.
    previous=array('I',(D[v//q][v % q] for v in range(q*q)))
    for length in range(3,k+1):
        old=previous
        previous=array('I',(old[v//q]*q+D[(v//q) % q][v % q] for v in range(q**length)))
    count=q**k; size=count*q
    with tempfile.TemporaryDirectory(prefix='kasa-latin-lift-',dir=pool.temp_dir) as d:
        root=Path(d)
        paths=[root/name for name in ('p.u32','c.u32','projection.u32','upper.u32','reserve.u32')]
        for path,data in zip(paths[:3],(P,C,previous)):
            with path.open('wb') as f: data.tofile(f)
        for path in paths[3:]:
            with path.open('wb') as f: f.truncate(size*4)
        parts=min(pool.workers*2,count)
        ranges=[(count*i//parts,count*(i+1)//parts) for i in range(parts)]
        tasks=[(*map(str,paths),q,R,D,start,stop) for start,stop in ranges]
        result=pool.map(_lift_range,tasks,progress=stage.pool_update if stage and stage.owner.enabled else None,task_units=[(stop-start)*q for start,stop in ranges])
        if result != ranges: raise RuntimeError('raw-lift worker range coverage differs')
        U=array('I'); V=array('I')
        with paths[3].open('rb') as f: U.fromfile(f,size)
        with paths[4].open('rb') as f: V.fromfile(f,size)
        return U,V,None
