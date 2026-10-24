"""Exact core solvers. Vertices must be 0,...,m-1; all k[v] are positive.

Tree optimization uses O(m) arithmetic operations. Optional counts enumerate
minimizing admissible CORE SETS, not all double Roman labellings.
"""
from math import inf


def adjacency(m, edges):
    adj=[[] for _ in range(m)]
    for u,v in edges:
        if u==v or not (0<=u<m and 0<=v<m):
            raise ValueError('Invalid edge')
        adj[u].append(v)
        adj[v].append(u)
    return adj


def check_set(adj, k, selected):
    if any(k[v]==1 and all(u in selected for u in adj[v]) for v in selected):
        return False
    if not selected:
        return True
    reached={next(iter(selected))}
    stack=list(reached)
    while stack:
        v=stack.pop()
        for u in adj[v]:
            if u in selected and u not in reached:
                reached.add(u)
                stack.append(u)
    return reached==selected


def tree_dp(adj, k, count=False, root=0):
    m=len(k)
    if m==0 or len(adj)!=m or min(k)<1:
        raise ValueError('Nonempty tree and positive multiplicities required')
    parent=[-2]*m
    parent[root]=-1
    order=[root]
    for v in order:
        for u in adj[v]:
            if u==parent[v]:
                continue
            if parent[u]!=-2:
                raise ValueError('Not a tree')
            parent[u]=v
            order.append(u)
    if len(order)!=m:
        raise ValueError('Disconnected tree')
    A=[inf]*m
    B=[inf]*m
    ca=[0]*m if count else None
    cb=[0]*m if count else None
    transitions=0
    for v in reversed(order):
        # d0: all processed children selected; d1: some child excluded.
        d0,d1=0,inf
        c0,c1=1,0
        for u in adj[v]:
            if u==parent[v]:
                continue
            transitions+=1
            nd0=d0+B[u]
            option1=d1+min(0,B[u])
            option2=d0  # exclude u for the first time
            nd1=min(option1,option2)
            if count:
                child_count=1 if B[u]>0 else cb[u] if B[u]<0 else 1+cb[u]
                nc0=c0*cb[u] if nd0<inf else 0
                nc1=0
                if option1==nd1 and nd1<inf:
                    nc1+=c1*child_count
                if option2==nd1 and nd1<inf:
                    nc1+=c0
                c0,c1=nc0,nc1
            d0,d1=nd0,nd1
        best=min(d0,d1)
        A[v]=k[v]-2+best
        B[v]=k[v]-2+(d1 if k[v]==1 else best)
        if count:
            ca[v]=(c0 if d0==best else 0)+(c1 if d1==best else 0)
            cb[v]=c1 if k[v]==1 else ca[v]
    optimum=0
    top=None
    number=1  # empty core set
    for v in order:
        value=B[v] if v==root else A[v]
        ways=(cb[v] if v==root else ca[v]) if count else 0
        if value<optimum:
            optimum,top,number=value,v,ways
        elif value==optimum and count:
            number+=ways
    # A canonical witness excludes zero-cost child states.
    selected=set()
    stack=[] if top is None else [(top,top!=root)]
    while stack:
        v,parent_outside=stack.pop()
        selected.add(v)
        children=[u for u in adj[v] if u!=parent[v]]
        chosen=[u for u in children if B[u]<0]
        if k[v]==1 and not parent_outside and len(chosen)==len(children):
            if not children:
                raise AssertionError('Impossible state chosen')
            chosen.remove(max(children,key=lambda u:B[u]))
        stack.extend((u,False) for u in chosen)
    assert check_set(adj,k,selected)
    assert sum(k[v]-2 for v in selected)==optimum
    return dict(correction=int(optimum),value=sum(k)+2*m+int(optimum),
                selected=selected,count=number if count else None,
                transitions=transitions,A=A,B=B,parent=parent)


def subsets(adj,k):
    m=len(k)
    masks=[sum(1<<u for u in ns) for ns in adj]
    best,count,witness=0,1,set()
    for mask in range(1,1<<m):
        vertices=[v for v in range(m) if mask>>v&1]
        if any(k[v]==1 and masks[v]&~mask==0 for v in vertices):
            continue
        reached=0
        front=mask&-mask
        while front:
            reached|=front
            neighbors=0
            while front:
                bit=front&-front
                front-=bit
                neighbors|=masks[bit.bit_length()-1]
            front=neighbors&mask&~reached
        if reached!=mask:
            continue
        value=sum(k[v]-2 for v in vertices)
        if value<best:
            best,count,witness=value,1,set(vertices)
        elif value==best:
            count+=1
    return dict(correction=best,value=sum(k)+2*m+best,count=count,selected=witness)


def leaf_response(adj,k,vertex):
    """rho after adding t>=0 private leaves at a vertex with k[v]>=2."""
    if k[vertex]<2:
        raise ValueError('The boundary family changes when a one-leaf core is modified')
    state=tree_dp(adj,k,root=vertex)
    outside=min([0]+[state['A'][u] for u in range(len(k)) if u!=vertex])
    inside=state['B'][vertex]
    return dict(inside=int(inside),outside=int(outside),
                threshold=max(0,int(outside-inside)),baseline=state['value'])


def rooted_milp(adj,k,time_limit=30):
    """One exact single-commodity-flow MILP for each possible selected root."""
    import numpy as np
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    m=len(k)
    arcs=[(u,v) for u,ns in enumerate(adj) for v in ns]
    dim=m+len(arcs)
    objective=np.zeros(dim)
    objective[:m]=np.asarray(k)-2
    best,witness=0,set()
    statuses=[]
    for root in range(m):
        lo=np.zeros(dim)
        hi=np.full(dim,max(0,m-1),dtype=float)
        hi[:m]=1
        lo[root]=1
        integral=np.zeros(dim)
        integral[:m]=1
        rows=[]
        lower=[]
        upper=[]
        def add(d,a,b):
            rows.append(d)
            lower.append(a)
            upper.append(b)
        for v in range(m):
            if k[v]==1:
                d={v:1}
                for u in adj[v]:
                    d[u]=1
                add(d,-np.inf,len(adj[v]))
        for j,(u,v) in enumerate(arcs):
            add({m+j:1,u:-(m-1)},-np.inf,0)
            add({m+j:1,v:-(m-1)},-np.inf,0)
        for v in range(m):
            d={v:-1} if v!=root else {u:-1 for u in range(m) if u!=root}
            for j,(a,b) in enumerate(arcs):
                sign=(int(a==v)-int(b==v)) if v==root else (int(b==v)-int(a==v))
                if sign:
                    d[m+j]=sign
            add(d,0,0)
        ri,ci,data=[],[],[]
        for i,row in enumerate(rows):
            for j,value in row.items():
                ri.append(i);ci.append(j);data.append(value)
        matrix=coo_matrix((data,(ri,ci)),shape=(len(rows),dim)).tocsc()
        result=milp(objective,integrality=integral,bounds=Bounds(lo,hi),
                    constraints=LinearConstraint(matrix,lower,upper),
                    options={'mip_rel_gap':0,'time_limit':time_limit})
        statuses.append(int(result.status))
        if result.status==2:
            continue
        if not result.success or result.mip_gap!=0:
            raise RuntimeError(str(result))
        selected={v for v in range(m) if result.x[v]>.5}
        assert check_set(adj,k,selected)
        value=sum(k[v]-2 for v in selected)
        assert abs(value-result.fun)<1e-6
        if value<best:
            best,witness=value,selected
    return dict(correction=best,value=sum(k)+2*m+best,selected=witness,statuses=statuses)
