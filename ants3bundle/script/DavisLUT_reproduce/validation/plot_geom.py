#!/usr/bin/env python3
# Compare the geometry-simulation reflected theta_out distribution (data/geom_<a>.txt, from
# validate_geom.js) against the LUT prediction, in linear and polar form. Writes to figures/.
import os, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__))
LUT=os.path.join(HERE,"..","luts","pooled_28um_fwd.lut")
FIG=os.path.join(HERE,"figures"); os.makedirs(FIG,exist_ok=True)
ANG=[10,45,80]
D=json.load(open(LUT)); B=D["Binning"]; nTi=B["ThetaIncBins"]; nTo=B["ThetaOutBins"]; nPh=B["PhiOutBins"]
def lut_marginal(a):
    x=a/(90.0/nTi)-0.5; k0=int(np.floor(x)); f=x-k0
    if k0<0: k0=0; f=0
    if k0>=nTi-1: k0=nTi-1; f=0
    k1=k0+1 if f>0 else k0
    def marg(k):
        h=np.array(D["ReflectedHist"][k],float).reshape(nTo,nPh).sum(axis=1); return h/h.sum() if h.sum()>0 else h
    m=(1-f)*marg(k0)+f*marg(k1); thc=(np.arange(nTo)+0.5)*90.0/nTo
    return thc, m/(90.0/nTo)
def meas(a):
    d=np.loadtxt(os.path.join(HERE,"data",f"geom_{a}.txt")); c=d[:,0]; n=d[:,1]/d[:,1].sum()
    return c, n/(c[1]-c[0])
# linear
fig,axes=plt.subplots(1,3,figsize=(14,4.2))
for i,a in enumerate(ANG):
    tl,ml=lut_marginal(a); tm,mm=meas(a)
    axes[i].plot(tl,ml,'-',color="tab:blue",lw=2,label="LUT prediction")
    axes[i].step(tm,mm,where="mid",color="tab:red",lw=1.2,label="geometry sim (detector)")
    axes[i].axvline(a,color="green",ls="--",lw=1,label=f"specular {a}°")
    axes[i].set_xlabel("θ_out from normal (deg)"); axes[i].set_ylabel("prob. density (/deg)")
    axes[i].set_title(f"{a}° incidence"); axes[i].grid(alpha=.3); axes[i].legend(fontsize=8)
fig.suptitle("Reflected θ_out: geometry simulation (source+monitor) vs LUT prediction",fontsize=13)
fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(os.path.join(FIG,"validation_geom_linear.png"),dpi=140)
# polar
fig=plt.figure(figsize=(13,4.6))
for i,a in enumerate(ANG):
    tl,ml=lut_marginal(a); tm,mm=meas(a)
    ax=fig.add_subplot(1,3,i+1,projection="polar")
    ax.set_thetamin(0); ax.set_thetamax(90); ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.plot(np.radians(tl),ml,'-',color="tab:blue",lw=2,label="LUT")
    ax.plot(np.radians(tm),mm,'o',color="tab:red",ms=2,label="geometry sim")
    ax.plot([np.radians(a)]*2,[0,max(ml.max(),mm.max())],color="green",ls="--",lw=1.2,label=f"specular {a}°")
    ax.set_title(f"{a}° incidence",fontsize=10,pad=12); ax.legend(fontsize=7,loc="upper right",bbox_to_anchor=(1.16,1.1))
fig.suptitle("Reflected θ_out in polar form: geometry simulation vs LUT",fontsize=12)
fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(FIG,"validation_geom_polar.png"),dpi=140)
print("wrote figures/validation_geom_linear.png and validation_geom_polar.png")
