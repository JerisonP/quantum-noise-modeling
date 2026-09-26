import numpy as np, time
from scipy.special import gamma
rng=np.random.default_rng(1)
def h_pulse(a,N):
    h=np.empty(N); h[0]=1
    for k in range(1,N): h[k]=h[k-1]*(a/2+k-1)/k
    return h
def fir(a,Qpsd,dt,N,M):
    Qd=Qpsd/dt**(1-a); h=h_pulse(a,N); nn=2*N
    H=np.fft.rfft(np.r_[h,np.zeros(N)])
    W=np.sqrt(Qd)*rng.standard_normal((M,N))
    return np.fft.irfft(H*np.fft.rfft(np.c_[W,np.zeros((M,N))],axis=1),nn,axis=1)[:,:N], Qd, h
def acf_demeaned(X,L):
    Y=X-X.mean(1,keepdims=True); N=X.shape[1]
    return np.array([np.mean(np.sum(Y[:,:N-m]*Y[:,m:],1)/(N-m)) for m in range(L+1)])
N=4001; dt=1e-3; Q=1e-3; M=400; L=int(0.3*N)
X,Qd,h=fir(1.0,Q,dt,N,M)
emp=acf_demeaned(X,L); rho=emp/emp[0]
lib=np.array([h[:N-m]@h[m:N] for m in range(L+1)]); lib/=lib[0]
# exact expectation of the estimator actually used
Hm=np.zeros((N,N))
for j in range(N): Hm[j:,j]=h[:N-j]
S=Qd*Hm@Hm.T
rm=S.mean(1,keepdims=True); S_dm=S-rm-rm.T+S.mean()
exact=np.array([np.mean(np.diagonal(S_dm,m)) for m in range(L+1)])
print("A) pink α=1 normalized-ACF RMSE vs LIBRARY theory: %.4f"%np.sqrt(np.mean((rho-lib)**2)))
print("   ... vs EXACT expectation of the demeaned estimator: %.4f"%np.sqrt(np.mean((rho-exact/exact[0])**2)))
# B) OU demean bias
tau=0.5; s=1; dto=0.01; To=10; n=int(To/dto)+1; Mo=400
phi=np.exp(-dto/tau); x=np.empty((Mo,n)); x[:,0]=rng.standard_normal(Mo)
for i in range(1,n): x[:,i]=phi*x[:,i-1]+np.sqrt(1-phi**2)*rng.standard_normal(Mo)
k=int(tau/dto); c=acf_demeaned(x,k)[k]
ckm=np.mean(np.sum(x[:,:n-k]*x[:,k:],1)/(n-k))
print("B) OU C(τc): demeaned %.3f | known-mean %.3f | theory %.3f | predicted bias 2τσ²/T=%.3f"%(c,ckm,np.exp(-1),2*tau/To))
# D) Kasdin gamma formula at large lag, α=0.5
a=0.5; m=np.arange(0,801)
with np.errstate(all='ignore'):
    R=(-1.0)**m*gamma(1-a)/(gamma(1+m-a/2)*gamma(1-m-a/2))
print("D) α=0.5 Gamma-formula ACF: first non-finite lag =", m[~np.isfinite(R)][0] if (~np.isfinite(R)).any() else None)
# E) WK rectangle-rule offset for OU
lags=np.arange(0,int(0.2*n)+1)*dto; C=np.exp(-lags/tau); f=np.array([0.0,5.0,20.0,49.0])
wk=[4*dto*np.sum(C*np.cos(2*np.pi*ff*lags)) for ff in f]; th=4*tau/(1+(2*np.pi*f*tau)**2)
print("E) WK (library, exact C input) / theory at f=",f,":",np.round(np.array(wk)/th,3))
