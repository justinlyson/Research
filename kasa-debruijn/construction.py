"""Revised-paper construction. Standard library only; internal infinity=q-1.

No search for Hamilton cycles: only paper-specified seeds, packet supplies,
and a deterministic spanning tree of their actual raw-cycle graph.
"""
import time
from array import array
from itertools import combinations


def require(ok, message):
    if not ok:
        raise ValueError(message)


def encode(w, q):
    v = 0
    for x in w:
        v = v*q+x
    return v


def word(v, q, k):
    w = [0]*k
    for i in range(k-1, -1, -1):
        v, w[i] = divmod(v, q)
    return tuple(w)


def shift(w, t, q):
    return tuple(q-1 if x == q-1 else (x+t) % (q-1) for x in w)


def rows(q, mode):
    n = q-1
    S = lambda x, t: n if x == n else (x+t) % n
    def rho(a, x):
        return x if a == n else a if x == n else n if x == a else (2*a-x) % n
    def W(x):
        m = n//2
        return 0 if x == n else n if x == m else n-1-x if x < m else n-x
    if mode == 'even':
        R = [[rho(a, S(x, 2)) for x in range(q)] for a in range(q)]
    elif mode in ('reverse', 'half'):
        t = n-2 if mode == 'reverse' else n//2
        R = [[S(x, t) if a == n else S(W(S(x, t-a)), a) for x in range(q)] for a in range(q)]
    else:
        R = [[(x+a) % q for x in range(q)] for a in range(q)]
    inv = [[0]*q for _ in range(q)]
    D = [[-1]*q for _ in range(q)]
    for a, row in enumerate(R):
        require(sorted(row) == list(range(q)), 'row is not a permutation')
        for x, y in enumerate(row):
            require(D[x][y] == -1, 'cover is not Latin')
            D[x][y] = a
            inv[a][y] = x
    return R, inv, D


def labels(P,stage=None):
    counted=0
    ids = array('I', [0xffffffff])*len(P)
    lengths = []
    for s in range(len(P)):
        if ids[s] != 0xffffffff:
            continue
        v, length = s, 0
        while ids[v] == 0xffffffff:
            ids[v] = len(lengths)
            v = P[v]
            length += 1
            if stage is not None and length & 0x3ffff==0:stage.update(counted+length)
        require(v == s, 'successor is not a permutation')
        lengths.append(length)
        counted+=length
    return ids, lengths


def complement(P, q, k):
    """Only used for memory-two seeds, never as a fallback constructor."""
    C = array('I')
    for v in range(len(P)):
        w = word(v, q, k)
        used = {shift((P[encode(shift(w, -t, q), q)] % q,), t, q)[0] for t in range(q-1)}
        require(len(used) == q-1, 'seed phases overlap')
        C.append((v % (q**(k-1)))*q+next(iter(set(range(q))-used)))
    require(len(set(C)) == len(C), 'seed complement is not a permutation')
    return C


def from_circular(s, q):
    a = list(map(int, s))
    require(len(a) == q*q, 'seed length')
    P = array('I', [0])*(q*q)
    seen = set()
    for i in range(len(a)):
        v = a[i]*q+a[(i+1) % len(a)]
        require(v not in seen, 'seed repeated vertex')
        seen.add(v)
        P[v] = a[(i+1) % len(a)]*q+a[(i+2) % len(a)]
    return P


def seed(q):
    n = q-1
    if q == 3:
        return from_circular('001122021', 3), None, {'branch': 'ternary-additive'}
    if q in (5, 7):
        P = from_circular({5:'0033421022011231413243044', 7:'0025114422346401530660316205541243261045213563365'}[q], q)
        return P, complement(P, q, 2), {'branch': 'finite-odd'}
    if q % 2 == 0:
        N = (q-2)//2
        if q in (4, 6):
            tau = {4:[0,2,1], 6:[0,2,3,1,4]}[q]
        else:
            tau = []
            for x in range(n):
                z = (x+N) % n-N
                tau.append((z if z >= 0 else -1 if z == -N else z-1) % n)
        def rho(a, x):
            return a if x == n else n if x == a else (2*a-x) % n
        C = array('I', (v*q+(u if v == n or u == n else (v+tau[(u-v) % n]) % n) for u in range(q) for v in range(q)))
        P = array('I', (v*q+(rho(1,u) if v == n else rho(0,C[u*q+v] % q)) for u in range(q) for v in range(q)))
        if q >= 8:
            ids, raw = labels(P)
            require(len(raw) == 2, 'negative-cycle seed must have two cycles')
            g = next((g for g in ([2,4] if N >= 4 else [2]) if ids[g*q] != ids[((1-g) % n)*q]), None)
            require(g is not None, 'paper offered seed pairs do not cross')
            x,y = g*q, ((1-g) % n)*q
            P[x],P[y] = P[y],P[x]
            for t in range(n):
                x,y = ((g+t) % n)*q+t, ((1-g+t) % n)*q+t
                C[x],C[y] = C[y],C[x]
        else:
            g = None
        require(len(labels(P)[1]) == 1, 'even seed is not Hamiltonian')
        require(C == complement(P,q,2), 'covariant seed residual correction differs')
        return P,C,{'branch':'finite-even' if q < 8 else 'uniform-even-negative-cycle', 'g':g}
    m,r = n//2,n//4
    W = lambda x: 0 if x == n else n if x == m else n-1-x if x < m else n-x
    F = lambda x: 0 if x == n else n if x == 0 else m-x if x < m else (m-1-x) % n
    P = array('I', (v*q+(F(u) if v == n else W(u)) for u in range(q) for v in range(q)))
    if q % 4 == 1:
        edits=[(n,m,0,m),(n-1,m,1,0),(m-1,m,m,1),(n,0,0,m),(m-1,0,m,0),(r-1,3*r,3*r,r+1),(3*r-1,3*r,r+1,3*r)]
    else:
        edits=[(n,m,0,m),(m-1,m,m,0),(n,0,0,m),(m-1,0,m,0),(3*r+2,r+1,r,3*r),(3*r+1,r+1,r+1,r),(r,r+1,3*r+1,r+1),(r+1,r+1,3*r,3*r+1),(r+1,3*r+2,3*r,3*r+1),(r,3*r+2,3*r+1,r+1),(3*r+1,3*r+2,r+1,r),(3*r+2,3*r+2,r,3*r),(0,m+1,n-1,n),(m,m+1,n,m-1),(m+1,m+1,m-1,n-1),(m+1,1,m-1,n-1),(m,1,n,m-1),(0,1,n-1,n)]
    for u,v,before,after in edits:
        require(P[u*q+v] % q == before, 'bounded recoloring old value differs')
        P[u*q+v] = v*q+after
    require(len(labels(P)[1]) == 1, 'bounded odd seed is not Hamiltonian')
    C = complement(P,q,2)
    return P,C,{'branch':'bounded-odd-seven-cell' if q % 4 == 1 else 'bounded-odd-eighteen-cell', 'edits':len(edits)}


def _raw_lift(P, C, q, k, R, D, pool=None,stage=None):
    if C is not None and pool is not None and pool.workers > 1 and q**(k+1) >= 2_000_000:
        from parallel_math import parallel_lift
        return parallel_lift(P,C,q,k,R,D,pool,stage)
    delta = array('I', range(q))
    for length in range(2,k+2):
        if length == 2:
            delta = array('I', (D[v//q][v % q] for v in range(q*q)))
        else:
            prev = delta
            delta = array('I', (prev[v//q]*q+D[(v//q) % q][v % q] for v in range(q**length)))
    mod = q**k
    def values(F,visible):
        for v,d in enumerate(delta):
            yield (v % mod)*q+R[F[d] % q][v % q]
            if visible and stage is not None and v & 0x3ffff==0:stage.update(v)
    U = array('I',values(P,C is None))
    V = None if C is None else array('I',values(C,True))
    return U,V,delta


def raw_lift(P,C,q,k,R,D,pool=None,progress=None):
    stage=progress.stage('Construction memory %d: raw Latin lift'%(k+1),q**(k+1),'vertices')if progress else None
    result=_raw_lift(P,C,q,k,R,D,pool,stage)
    if stage:stage.finish()
    return result


def partner(P,C,w,q):
    k = len(w)
    x = encode(w,q)
    for h in range(q):
        y = encode((h,)+w[1:],q)
        if x != y and P[x] == C[y] and P[y] == C[x]:
            return w,word(y,q,k)
    raise ValueError('required complete partner missing: '+str(w))


def geometric(R,a,s,k):
    w=[s]
    for _ in range(k-1):
        w.append(R[a][w[-1]])
    return tuple(w)


def packet(pair,q,R,inv):
    u,v=pair
    out=[]
    for y in range(q):
        tail=[y]
        for a in u[1:]:
            tail.append(R[a][tail[-1]])
        out.append((encode((inv[u[0]][y],*tail),q),encode((inv[v[0]][y],*tail),q)))
    return out


class Forest:
    def __init__(self,n): self.p=list(range(n)); self.count=n
    def find(self,x):
        while self.p[x] != x:
            self.p[x]=self.p[self.p[x]]; x=self.p[x]
        return x
    def union(self,a,b):
        a,b=self.find(a),self.find(b)
        if a == b: return False
        self.p[a]=b; self.count-=1; return True


def repair(P,C,q,k,edges,optional=None,progress=None):
    stage=progress.stage('Construction memory %d: raw-cycle labelling'%k,len(P),'vertices')if progress else None
    ids,lengths=labels(P,stage)
    if stage:stage.finish('%d actual raw cycles'%len(lengths))
    repair_start=time.perf_counter()
    require(len(lengths) <= q, 'raw lift has too many cycles')
    forest=Forest(len(lengths)); selected=[]
    for x,y in edges:
        if forest.union(ids[x],ids[y]): selected.append((x,y))
    geometric_components=forest.count
    extra=None
    if optional is not None:
        require(forest.count in (1,2), 'geometric graph has more than two components')
        if forest.count == 2:
            extra=next(((x,y) for x,y in optional if forest.find(ids[x]) != forest.find(ids[y])),None)
            require(extra is not None, 'no offered infinity occurrence crosses')
            edges=edges+[extra]
            x,y=extra
            forest.union(ids[x],ids[y]); selected.append(extra)
    require(forest.count == 1, 'paper packet graph disconnected')
    touched=set(); allphases=set()
    for x,y in edges:
        require(x != y and x//q**(k-1) != y//q**(k-1) and x % q**(k-1) == y % q**(k-1), 'packet suffix mismatch')
        require(P[x] == C[y] and P[y] == C[x], 'offered upper pair is not complete')
        require(x not in touched and y not in touched, 'packet endpoints overlap')
        touched.update((x,y))
    for t in range(q-1):
        part={encode(shift(word(v,q,k),t,q),q) for v in touched}
        require(not allphases.intersection(part), 'packet phase supports overlap')
        allphases.update(part)
    for x,y in selected: P[x],P[y]=P[y],P[x]
    for t in range(q-1):
        for x,y in selected:
            a=encode(shift(word(x,q,k),t,q),q); b=encode(shift(word(y,q,k),t,q),q)
            C[a],C[b]=C[b],C[a]
    if progress:
        progress.timings.append(dict(stage='Construction memory %d: repair'%k,group='construction',seconds=time.perf_counter()-repair_start))
        progress.log('Construction memory %d: repair complete; %d candidate pairs; %d selected exchanges; %d phases checked'%(k,len(edges),len(selected),q-1))
    return P,C,dict(memory=k,raw_cycle_lengths=lengths,offered_pairs=len(edges),tree_exchanges=len(selected),geometric_components=geometric_components,extra_infinity_pair=[word(x,q,k) for x in extra] if extra else None,
                    selected_pairs=[[word(x,q,k),word(y,q,k)] for x,y in selected],phase_supports_checked=q-1)


def check_geometric(P,C,q,k,R,a):
    for s in range(q):
        w=geometric(R,a,s,k); partner(P,C,w,q)
        nxt=encode(geometric(R,a,R[a][s],k),q)
        require(nxt in (P[encode(w,q)],C[encode(w,q)]), 'geometric mixed edge missing')


def even_renew(P,C,q,k,center,pool=None,progress=None):
    """Exactly three geometric packets and zero/one infinity occurrence."""
    R,inv,D=rows(q,'even')
    check_geometric(P,C,q,k,R,center)
    pairs=[partner(P,C,geometric(R,center,s,k),q) for s in (q-1,center,(center-2) % (q-1))]
    infinity=partner(P,C,(q-1,)*k,q)
    finite=partner(P,C,((center-1) % (q-1),)*k,q)
    edges=[edge for pair in pairs for edge in packet(pair,q,R,inv)]
    extra=packet(infinity,q,R,inv)
    U,V,_=raw_lift(P,C,q,k,R,D,pool,progress)
    old_anchor={v:(U[v],V[v]) for pair in packet(finite,q,R,inv) for v in pair}
    P,C,record=repair(U,V,q,k+1,edges,extra,progress)
    require(all((P[v],C[v]) == outputs for v,outputs in old_anchor.items()),'finite anchor raw outputs touched')
    center=(center-1) % (q-1)
    check_geometric(P,C,q,k+1,R,center)
    return P,C,center,record


def ternary_step(P,k,progress=None):
    q=3; R,_,D=rows(3,'additive')
    U,_,delta=raw_lift(P,None,q,k,R,D,progress=progress)
    primitive=bytearray(len(P)); v=0; s=0
    for _ in range(len(P)):
        primitive[v]=s; s=(s+v//(3**(k-1))) % 3; v=P[v]
    require(v == 0 and s == 0,'ternary primitive did not close')
    def neg(w): return tuple((-a) % 3 for a in w)
    def L(w):
        v=encode(w,3); return (w[0]-primitive[delta[v]]) % 3
    section=[((-u) % 3,)+tuple(i % 3 for i in range(k)) for u in range(3)]
    chosen=None
    for x,y in combinations(section,2):
        if L(x) == L(y) or L(neg(x)) == L(neg(y)): continue
        for j in range(3):
            a=tuple((z+j) % 3 for z in x); b=tuple((z+j) % 3 for z in y)
            if {L(a),L(b)} != {L(neg(a)),L(neg(b))}: chosen=(a,b); break
        if chosen: break
    require(chosen is not None,'ternary normalized ports failed')
    x,y=chosen; pairs=[(encode(x,3),encode(y,3)),(encode(neg(x),3),encode(neg(y),3))]
    require(len({v for pair in pairs for v in pair}) == 4,'ternary supports overlap')
    stage=progress.stage('Construction memory %d: raw-cycle labelling'%(k+1),len(U),'vertices')if progress else None
    ids,lengths=labels(U,stage); forest=Forest(len(lengths))
    if stage:stage.finish('%d actual raw cycles'%len(lengths))
    require(len(lengths) == 3,'ternary raw cycle count differs')
    for a,b in pairs:
        require(forest.union(ids[a],ids[b]),'ternary pair not a tree edge')
        U[a],U[b]=U[b],U[a]
    if progress:progress.log('Construction memory %d: ternary repair complete; 2 selected exchanges'%(k+1))
    return U,dict(memory=k+1,raw_cycle_lengths=lengths,tree_exchanges=2,normalized_ports=chosen)


def construct(q,k,max_states=100_000_000,pool=None,progress=None):
    require(type(q) is int and 3 <= q <= 255,'supported alphabet range is 3..255')
    require(type(k) is int and k >= 2,'revised prescribed-symmetry theorem requires k>=2; k=1 is unsupported')
    require(q**k <= min(max_states,0xffffffff),'state limit exceeded; increase --max-states deliberately')
    stage=progress.stage('Construction memory 2: seed',1,'seed')if progress else None
    P,C,info=seed(q); steps=[]
    if stage:stage.finish(info['branch'])
    if q == 3:
        for memory in range(2,k): P,record=ternary_step(P,memory,progress); steps.append(record)
        info['steps']=steps; return P,info
    R,inv,D=rows(q,'even' if q % 2 == 0 else 'reverse')
    center=0; h={5:3,7:4}.get(q,(q-1)//2+1); witness=(0,1)
    for memory in range(2,k):
        if progress:progress.log('Construction level %d/%d (q=%d)'%(memory+1,k,q))
        if memory == 2:
            if q in (4,6):
                fixed={4:[((0,1,3),(2,1,3))],6:[((5,0,3),(3,0,3))]}[q]
                edges=[tuple(encode(w,q) for w in pair) for pair in fixed]
            else:
                n=q-1; N=(q-2)//2; m=n//2
                if q % 2 == 0:
                    pairs=[((0,0),(n,0)),((2,0),(n-1,0)),((n-1,n-2),(1,n-2)),((0,2),((2-N) % n,2))]
                elif q in (5,7):
                    a,b={5:(2,0),7:(3,6)}[q]; pairs=[((a,a),(b,a))]
                else:
                    pairs=[((m,m),(n,m))]
                    pairs += [((0,m),(n-1,m))] if q % 4 == 1 else [((m-1,m+1),(m,m+1)),((m-1,1),(m,1))]
                for pair in pairs:
                    require(partner(P,C,pair[0],q)[1] == pair[1],'initializer lower pair differs')
                edges=[edge for pair in pairs for edge in packet(pair,q,R,inv)]
            U,V,_=raw_lift(P,C,q,memory,R,D,pool,progress)
            P,C,record=repair(U,V,q,memory+1,edges,progress=progress)
            if q % 2 == 0: check_geometric(P,C,q,3,R,0)
        elif q % 2 == 0:
            P,C,center,record=even_renew(P,C,q,memory,center,pool,progress)
        else:
            R,inv,D=rows(q,'half')
            x,y=witness
            witness_center=D[x][y]
            require(witness_center != q-1 and R[witness_center][y] == x,'new half-turn witness color differs')
            ws=[tuple((x if i % 2 == 0 else y) for i in range(memory)),tuple((y if i % 2 == 0 else x) for i in range(memory))]
            require(ws == [geometric(R,witness_center,x,memory),geometric(R,witness_center,y,memory)],'physical witnesses do not match recolored half-turn fibers')
            pairs=[partner(P,C,w,q) for w in ws]
            for i,w in enumerate(ws): require(encode(ws[1-i],q) in (P[encode(w,q)],C[encode(w,q)]),'odd alternating mixed edge missing')
            constant=partner(P,C,(h,)*memory,q)
            edges=[edge for pair in pairs for edge in packet(pair,q,R,inv)]
            U,V,_=raw_lift(P,C,q,memory,R,D,pool,progress)
            old_anchor={v:(U[v],V[v]) for pair in packet(constant,q,R,inv) for v in pair}
            P,C,record=repair(U,V,q,memory+1,edges,progress=progress)
            record['recomputed_half_turn_witness_center']=witness_center
            require(all((P[v],C[v]) == outputs for v,outputs in old_anchor.items()),'odd anchor raw outputs touched')
            check_geometric(P,C,q,memory+1,R,h)
            witness=(q-1,h); m=(q-1)//2; f=m//2 if m % 2 == 0 else (3*m-1)//2
            h=(h+f) % (q-1)
            partner(P,C,(h,)*(memory+1),q)
        steps.append(record)
    stage=progress.stage('Construction: final Hamilton check',len(P),'vertices')if progress else None
    require(len(labels(P,stage)[1]) == 1,'final reference is not Hamiltonian')
    if stage:stage.finish()
    info['steps']=steps
    return P,info
