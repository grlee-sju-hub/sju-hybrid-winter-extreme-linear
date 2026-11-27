#!/usr/bin/env python3
"""09-1 MLR selection (ERA5 pool, C3S grid): best-subset (k<=3) per predictand {EC_day,EW_day,CDHD}.
closed-form LOOCV (PRESS) -> LOO-ACC & LOO-MSSS ; VIF<=5 ; variants ALL(any lead) vs LEADONLY(lead>=1).
All series detrended + z-scored. Pool = predictor_accepted_era5.npz. Out: data/mlr_select_era5.xlsx.
NOTE: max-LOO-ACC over many subsets is mildly OPTIMISTIC (selection bias) -> read k=1->2->3 marginal gain."""
import numpy as np, pandas as pd, itertools, xarray as xr
from pathlib import Path
P=Path("."); D=P/"data"
KMAX=3; VIF_MAX=5.0; YR0,YR1=1979,2024; years=np.arange(YR0,YR1+1)
def detr1(x):
    x=np.asarray(x,float); n=len(x); m=np.isfinite(x)
    if m.sum()<4: return x-np.nanmean(x)
    c=np.polyfit(np.arange(n)[m],x[m],1); return x-(c[0]*np.arange(n)+c[1])
def zscore(x):
    s=np.nanstd(x); return (x-np.nanmean(x))/(s if s>0 else 1.0)
def pred_series(name):
    if name=="EC_day": d=xr.open_dataset(D/"ERA5_DJF_EC_index_1979-2024.nc"); y=d.year.values.astype(int); v=d.korea_nday.values.astype(float); d.close()
    elif name=="EW_day": d=xr.open_dataset(D/"ERA5_DJF_EW_index_1979-2024.nc"); y=d.year.values.astype(int); v=d.korea_nday.values.astype(float); d.close()
    elif name=="CDHD":
        d=xr.open_dataset(D/"ERA5_JA_HW_index_1979-2024.nc"); yh=d.year.values.astype(int); hwd=d.korea_nday.values.astype(float); d.close()
        sp=pd.read_csv(D/"SPEI"/"SPEIbase_v2.11_spei3_korea_JA_1901-2024.csv"); sp=sp[(sp.year>=1979)&(sp.year<=2024)]
        spr=sp.spei3_JA.values.astype(float); spei=(spr-np.nanmean(spr))/np.nanstd(spr,ddof=1)
        m=pd.DataFrame({"year":yh,"hwd":hwd}).merge(pd.DataFrame({"year":sp.year.values.astype(int),"spei":spei}),on="year")
        zH=np.log1p(m.hwd.values)/np.nanstd(np.log1p(m.hwd.values),ddof=1); y=m.year.values.astype(int); v=zH*(-m.spei.values)
    ymap={int(a):b for a,b in zip(y,v)}; return np.array([ymap.get(int(Y),np.nan) for Y in years])
def fit_loo(Pm,y):
    n=len(y); X=np.column_stack([np.ones(n),Pm]); XtXi=np.linalg.inv(X.T@X)
    beta=XtXi@(X.T@y); e=y-X@beta; h=np.einsum("ij,jk,ik->i",X,XtXi,X); eloo=e/(1-h)
    sse=np.sum(eloo**2); sst=np.sum((y-y.mean())**2); yloo=y-eloo
    acc=np.corrcoef(yloo,y)[0,1] if np.std(yloo)>1e-9 else np.nan
    return acc,1-sse/sst,1-np.sum(e**2)/sst,beta
d=np.load(D/"predictor_accepted_era5.npz",allow_pickle=True)
A_idx=d["index"].astype(str); A_var=d["var"].astype(str); A_lead=d["lead"].astype(int)
A_sign=d["sign"].astype(int); A_lat=d["lat_mean"]; A_lon=d["lon_mean"]; A_gid=d["gid"].astype(int); OBS=d["obs_ts"]
OBSz=np.stack([zscore(detr1(OBS[i])) for i in range(OBS.shape[0])])
def tag(i): return f"{A_var[i]}L{A_lead[i]}{'+' if A_sign[i]>0 else '-'}@({A_lat[i]:.0f},{A_lon[i]:.0f})"
rows=[]
for name in ["EC_day","EW_day","CDHD"]:
    y=detr1(pred_series(name))
    for variant,sel in [("ALL",np.ones(len(A_idx),bool)),("LEADONLY",A_lead>=1)]:
        pool=np.where((A_idx==name)&sel)[0]
        for k in range(1,KMAX+1):
            best=None
            for combo in itertools.combinations(pool,k):
                fin=np.isfinite(y).copy()
                for i in combo: fin=fin&np.isfinite(OBSz[i])
                yv=y[fin]; Pm=np.column_stack([OBSz[i][fin] for i in combo])
                if k>=2:
                    R=np.corrcoef(Pm.T)
                    try: mv=float(np.max(np.diag(np.linalg.inv(R))))
                    except np.linalg.LinAlgError: continue
                    if mv>VIF_MAX: continue
                else: mv=1.0
                acc,msss,r2,beta=fit_loo(Pm,yv)
                if (best is None) or (np.isfinite(acc) and acc>best[0]):
                    best=(acc,msss,r2,beta,mv,combo)
            if best is None: continue
            acc,msss,r2,beta,mv,combo=best
            rows.append(dict(index=name,variant=variant,k=k,loo_acc=round(float(acc),3),
                             loo_msss=round(float(msss),3),insample_r2=round(float(r2),3),max_vif=round(float(mv),2),
                             n_pool=len(pool),predictors=" + ".join(tag(i) for i in combo),
                             beta=",".join(f"{b:+.2f}" for b in beta[1:]),gids="|".join(str(A_gid[i]) for i in combo)))
res=pd.DataFrame(rows)
pd.set_option("display.width",240,"display.max_colwidth",80)
for name in ["EC_day","EW_day","CDHD"]:
    print(f"\n===== {name} =====",flush=True)
    print(res[res["index"]==name][["variant","k","loo_acc","loo_msss","max_vif","predictors"]].to_string(index=False),flush=True)
res.to_excel(D/"mlr_select_era5.xlsx",index=False)
print("\n[saved] mlr_select_era5.xlsx",flush=True); print("MLR_DONE",flush=True)
