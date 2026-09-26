import numpy as np
import os
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_acf_ou_gamma_wk.py')).read()
exec(_src.split('N=4001')[0])   # reuse h_pulse / fir helpers
def binned_slope(f,S,margin=0.3,nb=24):
    lf=np.log10(f); fr=lf[-1]-lf[0]; k=(lf>=lf[0]+margin*fr)&(lf<=lf[-1]-margin*fr)
    lf,S=lf[k],S[k]; e=np.linspace(lf.min(),lf.max(),nb+1); xs=[];ys=[]
    for b in range(nb):
        i=(lf>=e[b])&((lf<=e[b+1]) if b==nb-1 else (lf<e[b+1]))
        if i.any(): xs.append(lf[i].mean()); ys.append(np.log10(S[i].mean()))
    return -np.polyfit(xs,ys,1)[0]
N=4001; dt=1e-3; Q=1e-3
for a in [0.5,0.8,1.0,1.5]:
    X,Qd,h=fir(a,Q,dt,N,500)
    F=np.fft.rfft(X-X.mean(1,keepdims=True),axis=1)
    j=np.arange(1,N//2+1); f=j/(N*dt)
    P=(2*dt/N*np.abs(F[:,j])**2).mean(0)
    Hm=np.zeros((N,N))
    for c in range(N): Hm[c:,c]=h[:N-c]
    EP=2*dt/N*Qd*np.sum(np.abs(np.fft.rfft(Hm,axis=0)[j,:])**2,1)
    Sd=2*Qd*dt/np.abs(2*np.sin(np.pi*f*dt))**a
    print("α=%.1f  β̂=%.4f | target from asymptotic S(f): %.4f | target from EXACT finite-N E[periodogram]: %.4f"%(a,binned_slope(f,P),binned_slope(f,Sd),binned_slope(f,EP)))
