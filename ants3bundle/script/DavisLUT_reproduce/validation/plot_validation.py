#!/usr/bin/env python3
# Rule-level validation, polar view: measured (real runtime sampling, from validate_lut) vs
# stored-LUT reflected theta_out distribution. Reads data/valid_<a>.txt, writes
# figures/validation_polar.png. The 2D (theta_out x phi_out) lego overlay is plot_lego.py.
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__))
FIG=os.path.join(HERE,"figures"); os.makedirs(FIG,exist_ok=True)
ANG=[10,45,80]
def load(a):
    fn=os.path.join(HERE,"data",f"valid_{a:02d}.txt"); hdr=open(fn).readline()
    meta=dict(t.split("=") for t in hdr.replace("# ","").split() if "=" in t)
    nTo=int(meta["nThetaOut"]); nPh=int(meta["nPhiOut"])
    d=np.loadtxt(fn); return d[:,2].reshape(nTo,nPh), d[:,3].reshape(nTo,nPh), nTo, nPh, float(meta["R_meas"]), float(meta["R_lut"])
fig=plt.figure(figsize=(13,4.6))
for i,a in enumerate(ANG):
    meas,pred,nTo,nPh,Rm,Rl=load(a)
    thc=(np.arange(nTo)+0.5)*90.0/nTo
    m=meas.sum(axis=1); p=pred.sum(axis=1); tvd=0.5*np.abs(m-p).sum()
    ax=fig.add_subplot(1,3,i+1,projection="polar")
    ax.set_thetamin(0); ax.set_thetamax(90); ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.plot(np.radians(thc),p,'-',color="tab:blue",lw=2,label="LUT prediction")
    ax.plot(np.radians(thc),m,'o',color="#ff1111",ms=3,label="measured")
    ax.plot([np.radians(a)]*2,[0,max(m.max(),p.max())],color="green",ls="--",lw=1.2,label=f"specular {a}°")
    ax.set_title(f"{a}° incidence  (R={Rm:.3f}/{Rl:.3f}, TVD={tvd:.3f})",fontsize=9,pad=12)
    ax.legend(fontsize=7,loc="upper right",bbox_to_anchor=(1.15,1.1))
fig.suptitle("Rule-level validation: measured vs stored LUT (reflected θ_out, polar)",fontsize=12)
fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(FIG,"validation_polar.png"),dpi=140)
print("wrote figures/validation_polar.png")
