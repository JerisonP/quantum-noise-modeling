# Statistical tests of the generators, run over 200 seeds each: count how often
# a 5-SE check fails (should be ~0) and report the worst |z|.
import numpy as np
def h_pulse(a,N):
    h=np.empty(N); h[0]=1
    for k in range(1,N): h[k]=h[k-1]*(a/2+k-1)/k
    return h
def frac_gen(a,Q,dt,N,M,rng):
    Qd=Q*dt**(a-1); h=h_pulse(a,N); W=np.sqrt(Qd)*rng.standard_normal((N,M))
    nf=2*N; H=np.fft.rfft(h,nf); return np.fft.irfft(H[:,None]*np.fft.rfft(W,nf,axis=0),nf,axis=0)[:N], Qd*np.sum(h**2)
def ou_gen(mu,s,tau,dt,n,M,rng):
    phi=np.exp(-dt/tau); sd=s*np.sqrt(-np.expm1(-2*dt/tau)); X=np.empty((n,M)); x=s*rng.standard_normal(M); X[0]=mu+x
    for k in range(1,n): x=phi*x+sd*rng.standard_normal(M); X[k]=mu+x
    return X
zs={}
for seed in range(200):
    rng=np.random.default_rng(seed)
    X,v=frac_gen(1.0,1e-3,1e-3,512,4000,rng); M=4000
    zs.setdefault('frac var',[]).append((X[-1].var(ddof=1)-v)/(v*np.sqrt(2/(M-1))))
    mu,s,tau,dt,M=0.7,2.0,0.5,0.01,4000
    X=ou_gen(mu,s,tau,dt,301,M,rng)
    zs.setdefault('ou var',[]).append((X[150].var(ddof=1)-s*s)/(s*s*np.sqrt(2/(M-1))))
    p=(X[50]-mu)*(X[100]-mu); zs.setdefault('ou C(τc)',[]).append((p.mean()-s*s/np.e)/(p.std(ddof=1)/np.sqrt(M)))
    y0,y1=X[99]-mu,X[100]-mu; ph=np.sum(y0*y1)/np.sum(y0**2)
    zs.setdefault('ou ϕ',[]).append((ph-np.exp(-dt/tau))/np.sqrt((1-np.exp(-2*dt/tau))/M))
for k,z in zs.items():
    z=np.array(z); print(f"{k:10s} mean z={z.mean():+.2f} sd z={z.std():.2f} max|z|={np.abs(z).max():.2f} fails(>5)={np.sum(np.abs(z)>5)}")
