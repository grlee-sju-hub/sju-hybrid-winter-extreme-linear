#!/usr/bin/env python3
"""Direct hybrid verify of the FIXED 09-1 ALL k3 models (no re-selection): L0 predictor -> C3S MME at predict time.
Two training protocols, both eval on 1993-2016 (C3S overlap):
  A) 1993-train : fit & LOO within 1993-2016 (detrend over 1993-2016). [apples-to-apples, keep]
  B) 1979-train : fit beta on 1979-2024 observed (LOO leaves out the eval year -> 45yr fit); predict 1993-2016.
     L0 C3S calibrated to observed-L0 distribution over the 1993-2016 overlap (mean+std match; period scales differ).
Lead terms always observed; coeffs always fit on OBSERVED. Out: console + data/mlr_pp_fixed_era5.csv."""
import numpy as np, pandas as pd, xarray as xr
from pathlib import Path
D=Path("data")
def detr(x):
    x=np.asarray(x,float); n=len(x); t=np.arange(n)-(n-1)/2.; a=x-x.mean(); s=(t*t).sum(); return a-((t*a).sum()/s)*t
def corr(a,b): return float(np.corrcoef(a,b)[0,1])
def pploo(Xf,Xp,y):
    XtXi=np.linalg.inv(Xf.T@Xf); beta=XtXi@(Xf.T@y); e=y-Xf@beta
    h=np.einsum("ij,jk,ik->i",Xf,XtXi,Xf); cross=np.einsum("ij,jk,ik->i",Xp,XtXi,Xf)
    yloo=Xp@beta-cross*e/(1-h); return corr(yloo,y)
def pred_series(name):
    if name=="EC_day": d=xr.open_dataset(D/"ERA5_DJF_EC_index_1979-2024.nc"); y=d.year.values.astype(int); v=d.korea_nday.values.astype(float); d.close()
    elif name=="EW_day": d=xr.open_dataset(D/"ERA5_DJF_EW_index_1979-2024.nc"); y=d.year.values.astype(int); v=d.korea_nday.values.astype(float); d.close()
    elif name=="CDHD":
        d=xr.open_dataset(D/"ERA5_JA_CDHD_index_1979-2024.nc"); y=d.year.values.astype(int); v=d.korea_nday.values.astype(float); d.close()
    return {int(a):float(b) for a,b in zip(y,v)}
d=np.load(D/"predictor_accepted_era5.npz",allow_pickle=True)
G=d["gid"].astype(int); OBS=d["obs_ts"]; LEAD=d["lead"].astype(int); OBSy={int(y):i for i,y in enumerate(d["years"].astype(int))}
gi={int(g):i for i,g in enumerate(G)}
mme=pd.read_csv(D/"c3s_mme_L0_djf.csv"); MME={}
for _,r in mme.iterrows(): MME.setdefault(int(r.gid),{})[int(r.year)]=float(r.value)
MODELS={"EC_day":[10,191,132],"EW_day":[218,313,265]}
Y46=np.arange(1979,2025); E93=np.arange(1993,2017)
rows=[]
for nm,gids in MODELS.items():
    ym=pred_series(nm); L0=[g for g in gids if LEAD[gi[g]]==0]
    # ---- A) 1993-train ----
    i93=[OBSy[Y] for Y in E93]; yA=detr(np.array([ym[Y] for Y in E93]))
    Xo=np.column_stack([np.ones(len(E93))]+[detr(OBS[gi[g]][i93]) for g in gids])
    Xh=np.column_stack([np.ones(len(E93))]+[(detr(np.array([MME[g][Y] for Y in E93])) if g in L0 else detr(OBS[gi[g]][i93])) for g in gids])
    accA_obs=pploo(Xo,Xo,yA); accA_hyb=pploo(Xo,Xh,yA)
    # ---- B) 1979-train, predict 1993-2016 (explicit LOO) ----
    y46=detr(np.array([ym[Y] for Y in Y46])); Xf46={g:detr(OBS[gi[g]][[OBSy[Y] for Y in Y46]]) for g in gids}
    pos={int(Y):k for k,Y in enumerate(Y46)}
    cal={}
    for g in L0:                                   # calibrate C3S -> obs-L0 distribution over overlap
        oL0=np.array([Xf46[g][pos[Y]] for Y in E93]); c=np.array([MME[g][Y] for Y in E93])
        cc=(c-c.mean())/c.std()*oL0.std()+oL0.mean(); cal[g]={int(Y):cc[k] for k,Y in enumerate(E93)}
    yo=[]; yh=[]; yt=[]
    for Y in E93:
        tr=Y46!=Y; Xtr=np.column_stack([np.ones(tr.sum())]+[Xf46[g][tr] for g in gids]); beta=np.linalg.lstsq(Xtr,y46[tr],rcond=None)[0]
        iY=pos[Y]; xo=np.array([1.0]+[Xf46[g][iY] for g in gids]); xh=np.array([1.0]+[(cal[g][Y] if g in L0 else Xf46[g][iY]) for g in gids])
        yo.append(xo@beta); yh.append(xh@beta); yt.append(y46[iY])
    accB_obs=corr(yo,yt); accB_hyb=corr(yh,yt)
    rows.append(dict(index=nm,nL0=len(L0),
        A_obs93=round(accA_obs,3),A_hyb93_C3S=round(accA_hyb,3),
        B_obs_tr79=round(accB_obs,3),B_hyb_tr79_C3S=round(accB_hyb,3)))
res=pd.DataFrame(rows); print(res.to_string(index=False),flush=True); res.to_csv(D/"mlr_pp_fixed_era5.csv",index=False)
print("[saved] mlr_pp_fixed_era5.csv  (eval 1993-2016; A=train93, B=train79)",flush=True)
