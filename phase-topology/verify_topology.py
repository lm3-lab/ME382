import math
R = 8.314

def envelope(T, p, N=401):
    TmA, TmB, dHA, dHB, Oml, Oms = p
    dx = 1.0/(N-1); pts=[]
    for i in range(N):
        x = i*dx
        s = 0.0 if (x<=0.0 or x>=1.0) else x*math.log(x)+(1-x)*math.log(1-x)
        Gl = (1-x)*dHA*(1-T/TmA) + x*dHB*(1-T/TmB) + Oml*x*(1-x) + R*T*s
        Gs = Oms*x*(1-x) + R*T*s
        pts.append((x,Gl,'L') if Gl < Gs else (x,Gs,'S'))
    return pts, dx

def lower_hull(pts):
    h=[]
    for q in pts:
        while len(h)>=2:
            (x0,y0,_),(x1,y1,_) = h[-2],h[-1]
            if (x1-x0)*(q[1]-y0)-(y1-y0)*(q[0]-x0) <= 0: h.pop()
            else: break
        h.append(q)
    return h

def state(T, p, N=401):
    pts,dx = envelope(T,p,N); h = lower_hull(pts)
    tl = [(a,b) for a,b in zip(h,h[1:]) if b[0]-a[0] > 1.5*dx]
    ss = [t for t in tl if t[0][2]=='S' and t[1][2]=='S']
    lc = sorted([a[0] for a,b in tl if a[2]=='L' and b[2]=='S'] +
                [b[0] for a,b in tl if a[2]=='S' and b[2]=='L'])
    return ss, lc

def find_invariant(p, THI, TLO, steps=900, N=401):
    """Cool until the first S|S tie-line. The liquid contact is read one step
    ABOVE it, where it is still well conditioned. No liquid contact there means
    the solid gap opened below the solidus -> lens, not an invariant."""
    prev_lc = None
    for k in range(steps+1):
        T = THI + (TLO-THI)*k/steps
        ss, lc = state(T,p,N)
        if ss:
            if not prev_lc: return None                  # lens
            xa, xb = ss[0][0][0], ss[0][1][0]
            return T, xa, xb, sum(prev_lc)/len(prev_lc)
        prev_lc = lc
    return None

def run(p, THI, TLO, tol=1.5/400):
    r = find_invariant(p,THI,TLO)
    if r is None: return ('lens',None,None,None,None)
    T,xa,xb,xL = r
    if abs(xL-xa) <= tol or abs(xL-xb) <= tol: k='xing'
    elif xa < xL < xb:                          k='eu'
    else:                                       k='PERI'
    return (k,T,xa,xb,xL)

def line(lbl,p,THI,TLO):
    k,T,xa,xb,xL = run(p,THI,TLO)
    if k=='lens': return f"{lbl}  lens"
    return f"{lbl}  {k:>5}  T={T:6.1f}  x_a={xa:.3f}  x_L={xL:.3f}  x_b={xb:.3f}"

print("=== Branch 1: TmA=700 TmB=900 dHA=dHB=12 ===")
for O in range(8,17):
    print(line(f"Om_s={O:2d}", (700.,900.,12000.,12000.,0.,O*1000.), 1000., 200.))

print("\n=== Branch 2: TmA=700 TmB=1400 dHA=12 dHB=40 ===")
for O in range(10,26):
    print(line(f"Om_s={O:2d}", (700.,1400.,12000.,40000.,0.,O*1000.), 1500., 200.))

print("\n=== Negative result: symmetric dHm=12 kJ/mol, vary TmB ===")
print("TmB  " + " ".join(f"{O:>5d}" for O in range(12,38,2)))
for TmB in (900.,1100.,1400.,1800.,2400.):
    row=[]
    for O in range(12,38,2):
        k,_,_,_,_ = run((700.,TmB,12000.,12000.,0.,O*1000.), max(1.4*TmB,1200.), 200.)
        row.append(k)
    print(f"{TmB:4.0f} " + " ".join(f"{c:>5}" for c in row))

print("\n=== dS_m asymmetry ===")
for nm,(dH,Tm) in (('A',(12000.,700.)),('B sym',(12000.,900.)),('B branch2',(40000.,1400.))):
    print(f"  {nm:10s} dS_m = {dH/Tm:5.1f} J/mol/K")
