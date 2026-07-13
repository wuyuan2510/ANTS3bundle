#!/usr/bin/env python3
# 3D lego/surface comparison of the reflected angular distribution (theta_out x phi_out):
# LUT prediction (blue surface) and measured/validation (red wireframe) overlaid in the SAME
# axes, one panel per incidence angle, log z-scale for readability. Reads data/valid_<a>.txt,
# writes figures/validation_heatmap.png.
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: registers 3d projection
HERE=os.path.dirname(os.path.abspath(__file__))
FIG=os.path.join(HERE,"figures"); os.makedirs(FIG,exist_ok=True)
ANG=[10,45,80]
def load(a):
    fn=os.path.join(HERE,"data",f"valid_{a:02d}.txt"); hdr=open(fn).readline()
    meta=dict(t.split("=") for t in hdr.replace("# ","").split() if "=" in t)
    nTo=int(meta["nThetaOut"]); nPh=int(meta["nPhiOut"])
    d=np.loadtxt(fn); return d[:,2].reshape(nTo,nPh), d[:,3].reshape(nTo,nPh), nTo, nPh
fig=plt.figure(figsize=(16,5.4))
for i,a in enumerate(ANG):
    meas,pred,nTo,nPh=load(a)
    phi=(np.arange(nPh)+0.5)*360.0/nPh; tho=(np.arange(nTo)+0.5)*90.0/nTo
    PHI,THO=np.meshgrid(phi,tho); FLOOR=1e-2
    Zl=np.log10(np.maximum(pred*100,FLOOR)); Zm=np.log10(np.maximum(meas*100,FLOOR))
    ax=fig.add_subplot(1,3,i+1,projection="3d")
    ax.plot_surface(PHI,THO,Zl,cmap="Blues",alpha=0.65,linewidth=0,antialiased=True)
    ax.plot_wireframe(PHI,THO,Zm,color="#ff1111",linewidth=0.7,rstride=1,cstride=2)
    ax.set_zlim(np.log10(FLOOR), max(Zl.max(),Zm.max())+0.1)
    ax.set_xlabel("φ_out from incidence plane (deg)",fontsize=8,labelpad=2)
    ax.set_ylabel("θ_out from normal (deg)",fontsize=8,labelpad=2)
    ax.set_zlabel("log10(probability/bin, %)",fontsize=8,labelpad=2)
    ax.set_title(f"{a}° incidence",fontsize=11); ax.view_init(elev=28,azim=-58); ax.tick_params(labelsize=7)
handles=[Line2D([0],[0],color="tab:blue",lw=6,alpha=0.65,label="LUT prediction (surface)"),
         Line2D([0],[0],color="#ff1111",lw=2.0,label="validation / measured (wireframe)")]
fig.suptitle("Reflected angular distribution — LUT vs validation overlaid (3D lego, log z)",y=0.98,fontsize=13)
fig.legend(handles=handles,loc="lower center",ncol=2,fontsize=10,frameon=True,bbox_to_anchor=(0.5,0.005))
fig.tight_layout(rect=(0,0.07,1,0.95)); fig.savefig(os.path.join(FIG,"validation_heatmap.png"),dpi=140)
print("wrote figures/validation_heatmap.png")
