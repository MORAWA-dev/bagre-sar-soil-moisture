"""Monthly SAR/optical predictors and SMAP-scale retrieval. CPU only."""
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.transform import Affine
from rasterio.warp import reproject, Resampling, transform_bounds
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import HistGradientBoostingRegressor

REGION=(-1.20,10.85,-0.05,12.05)
CORE=(-1.05,11.00,-0.20,11.90)
PILOT=(-0.59,11.41,-0.47,11.52)
CRS='EPSG:32630'
MONTHS=pd.date_range('2024-01-01','2025-12-01',freq='MS')
S1=['vv_power','vh_power','angle','sar_count']
S2=['ndvi','optical_count']
FEATURES=['vv_power','vh_power','angle','ndvi']
NODATA=-9999.0

def power(db): return 10.0**(np.asarray(db)/10.0)
def db(p): return 10*np.log10(np.maximum(p,1e-12))
def grid(bounds,crs,scale):
    l,b,r,t=transform_bounds('EPSG:4326',crs,*bounds,densify_pts=21)
    l,b=np.floor(np.array([l,b])/scale)*scale
    r,t=np.ceil(np.array([r,t])/scale)*scale
    return dict(crs=crs,transform=[scale,0,float(l),0,-scale,float(t)],width=int(round((r-l)/scale)),height=int(round((t-b)/scale)))

def native_grid(proj,bounds=REGION):
    """Subset the EE ingested SMAP grid without pretending it is original EASE-Grid."""
    crs=proj['crs']; a=Affine(*proj['transform'])
    assert abs(a.b)<1e-12 and abs(a.d)<1e-12 and a.a>0 and a.e<0, proj
    l,b,r,t=transform_bounds('EPSG:4326',crs,*bounds,densify_pts=21)
    c0,r0=(~a)*(l,t); c1,r1=(~a)*(r,b)
    c0,r0=int(np.floor(c0)),int(np.floor(r0)); c1,r1=int(np.ceil(c1)),int(np.ceil(r1))
    sub=a*Affine.translation(c0,r0)
    return dict(crs=crs,transform=list(sub)[:6],width=c1-c0,height=r1-r0)

def shape(g): return (g['height'],g['width'])

def write_tif(path,arr,g,bands):
    arr=np.asarray(arr,dtype=np.float32)
    if arr.ndim==2: arr=arr[None,...]
    assert arr.shape==(len(bands),*shape(g)), (arr.shape,g)
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix('.partial.tif')
    with rio.open(tmp,'w',driver='GTiff',width=g['width'],height=g['height'],count=len(bands),dtype='float32',crs=g['crs'],transform=Affine(*g['transform']),nodata=NODATA,compress='deflate',photometric='MINISBLACK') as dst:
        dst.write(np.where(np.isfinite(arr),arr,NODATA)); dst.descriptions=tuple(bands)
    tmp.replace(p)

def read_tif(path,g,bands):
    with rio.open(path) as src:
        assert src.shape==shape(g) and src.count==len(bands), (path,src.shape,src.count)
        assert src.crs==rio.crs.CRS.from_user_input(g['crs']), (path,src.crs)
        assert src.transform.almost_equals(Affine(*g['transform']),precision=1e-6), (path,src.transform)
        assert tuple(src.descriptions)==tuple(bands), (path,src.descriptions)
        return src.read(masked=True).filled(np.nan).astype('float32')

def download(image,path,g,bands,force=False):
    """Download on the exact grid, subdividing requests on EE reprojection limits."""
    import requests
    from rasterio.io import MemoryFile
    path=Path(path)
    if path.exists() and not force:
        read_tif(path,g,bands); print('Cached',path.name); return
    prepared=image.select(bands).toFloat().unmask(NODATA,False)
    arr=np.full((len(bands),*shape(g)),np.nan,dtype=np.float32)
    base=Affine(*g['transform'])
    def fetch(row,col,height,width):
        transform=base*Affine.translation(col,row)
        try:
            url=prepared.getDownloadURL({'crs':g['crs'],'crs_transform':list(transform)[:6],'dimensions':[width,height],'format':'GEO_TIFF','filePerBand':False})
            with requests.get(url,timeout=240) as response:
                if not response.ok:
                    raise RuntimeError(f'HTTP {response.status_code}: {response.text[:2000]}')
                payload=response.content
        except Exception as exc:
            message=str(exc).lower()
            if (any(marker in message for marker in ('too large','memory limit','memory capacity exceeded','timed out')) and max(height,width)>32):
                print(f'{path.name}: subdividing {width}x{height} request',flush=True)
                if width>=height:
                    mid=width//2; fetch(row,col,height,mid); fetch(row,col+mid,height,width-mid)
                else:
                    mid=height//2; fetch(row,col,mid,width); fetch(row+mid,col,height-mid,width)
                return
            raise
        with MemoryFile(payload) as mem, mem.open() as src:
            assert src.shape==(height,width) and src.count==len(bands), (src.shape,src.count,height,width)
            assert src.crs==rio.crs.CRS.from_user_input(g['crs'])
            assert src.transform.almost_equals(transform,precision=1e-6), (src.transform,transform)
            tile=src.read().astype('float32'); tile[tile==NODATA]=np.nan
            arr[:,row:row+height,col:col+width]=tile
    if max(shape(g))>200:
        for row in range(0,g['height'],200):
            for col in range(0,g['width'],200):
                fetch(row,col,min(200,g['height']-row),min(200,g['width']-col))
    else:
        fetch(0,0,*shape(g))
    write_tif(path,arr,g,bands); print(path.name,shape(g))

def aggregate(a,src,dst,method=Resampling.average):
    a=np.asarray(a,dtype=np.float32); out=np.full(shape(dst),NODATA,dtype=np.float32)
    reproject(np.where(np.isfinite(a),a,NODATA),out,src_transform=Affine(*src['transform']),src_crs=src['crs'],src_nodata=NODATA,dst_transform=Affine(*dst['transform']),dst_crs=dst['crs'],dst_nodata=NODATA,resampling=method)
    out[out==NODATA]=np.nan; return out

def centers(g):
    rr,cc=np.indices(shape(g)); a=Affine(*g['transform'])
    x,y=a*(cc+.5,rr+.5)
    return x,y

def core_mask(g):
    from pyproj import Transformer
    x,y=centers(g); lon,lat=Transformer.from_crs(g['crs'],'EPSG:4326',always_xy=True).transform(x,y)
    return (lon>=CORE[0])&(lon<=CORE[2])&(lat>=CORE[1])&(lat<=CORE[3])

def summarize_month(sar,optical,land,smap,fine,coarse,month):
    """Common-support coarse summaries. Coverage denominator includes water and missing areas."""
    valid=(land[0]>=.8)&np.isfinite(sar[:3]).all(0)&np.isfinite(optical[0])&(sar[3]>=1)&(optical[1]>=1)
    coverage=aggregate(valid.astype('float32'),fine,coarse)
    land_fraction=aggregate(land[0],fine,coarse)
    crop_fraction=aggregate(land[1],fine,coarse)
    bands=[aggregate(np.where(valid,v,np.nan),fine,coarse) for v in [*sar[:3],optical[0]]]
    rr,cc=np.indices(shape(coarse)); x,y=centers(coarse)
    data=dict(month=str(month)[:7],row=rr.ravel(),col=cc.ravel(),x=x.ravel(),y=y.ravel(),coverage=coverage.ravel(),land_fraction=land_fraction.ravel(),crop_fraction=crop_fraction.ravel(),sm=smap[0].ravel(),smap_days=smap[1].ravel(),core=core_mask(coarse).ravel())
    data.update({k:v.ravel() for k,v in zip(FEATURES,bands)})
    table=pd.DataFrame(data)
    table['predictable']=table.core&(table.coverage>=.75)&(table.land_fraction>=.8)&np.isfinite(table[FEATURES]).all(axis=1)
    table['eligible']=table.predictable&(table.smap_days>=4)&table.sm.between(0,.65)
    return table,valid

def assign_splits(df):
    """Spatial holdout is a contiguous northeast quadrant; temporal is entire 2025."""
    df=df.copy(); core=df[df.core]; assert len(core)>0
    row_cut=core.row.median(); col_cut=core.col.median()
    spatial=(df.row<=row_cut)&(df.col>=col_cut)
    year=pd.to_datetime(df.month).dt.year; month=pd.to_datetime(df.month).dt.month
    df['split']=np.select([spatial&(year==2025),spatial&(year==2024),~spatial&(year==2025),~spatial&(year==2024)&month.isin([4,8,12])],['joint_test','spatial_test','temporal_test','validation'],default='train')
    return df

def add_season(df):
    df=df.copy(); phase=2*np.pi*(pd.to_datetime(df.month).dt.month-1)/12
    df['sin_month']=np.sin(phase); df['cos_month']=np.cos(phase); return df

def candidates():
    season=['sin_month','cos_month']
    return {
      'seasonal_ridge':(season,make_pipeline(StandardScaler(),Ridge(alpha=5))),
      'sar_ridge':(FEATURES[:3],make_pipeline(StandardScaler(),Ridge(alpha=5))),
      'optical_ridge':(['ndvi'],make_pipeline(StandardScaler(),Ridge(alpha=5))),
      'fusion_ridge':(FEATURES,make_pipeline(StandardScaler(),Ridge(alpha=5))),
      'fusion_hgb':(FEATURES,HistGradientBoostingRegressor(max_iter=100,max_leaf_nodes=15,min_samples_leaf=20,l2_regularization=1,random_state=42))}

def metrics(y,p):
    y,p=np.asarray(y),np.asarray(p); n=len(y)
    if not n: return dict(n=0,bias=np.nan,mae=np.nan,rmse=np.nan,ubrmse=np.nan,r=np.nan,r2=np.nan)
    e=p-y; b=e.mean(); return dict(n=n,bias=b,mae=mean_absolute_error(y,p),rmse=np.sqrt(mean_squared_error(y,p)),ubrmse=np.sqrt(np.mean((e-b)**2)),r=np.corrcoef(y,p)[0,1] if n>2 and np.std(y)>0 and np.std(p)>0 else np.nan,r2=r2_score(y,p) if n>1 and np.std(y)>0 else np.nan)

def fit_models(df):
    df=add_season(df[df.eligible]); train=df[df.split=='train']; val=df[df.split=='validation']
    assert len(train)>=100 and len(val)>=25, 'Insufficient training/validation support; inspect audit, do not silently relax masks.'
    models={}; rows=[]; preds=[]
    for name,(cols,model) in candidates().items():
        model.fit(train[cols],train.sm); models[name]=(cols,model)
    for name,(cols,model) in {'train_mean':([],None),**models}.items():
        for split in ['validation','spatial_test','temporal_test','joint_test']:
            d=df[df.split==split].copy(); p=(np.full(len(d),train.sm.mean()) if model is None else model.predict(d[cols])) if len(d) else np.array([])
            rows.append(dict(model=name,split=split,**metrics(d.sm,p)))
            d['model']=name; d['prediction']=p; preds.append(d)
    scores=pd.DataFrame(rows); sel=scores[(scores.split=='validation')&scores.model.isin(models)].sort_values('rmse').iloc[0].model
    return models,scores,pd.concat(preds,ignore_index=True),sel

def fine_contrast(model,sar,optical,valid,fine,coarse,eligible):
    """Linear ridge score centered at each 250m pixel's coarse-grid mean; NOT soil moisture."""
    x=np.stack([*sar[:3],optical[0]],axis=-1); raw=np.full(valid.shape,np.nan,dtype=np.float32)
    if valid.any(): raw[valid]=model.predict(pd.DataFrame(x[valid],columns=FEATURES))
    mean=aggregate(raw,fine,coarse); centered=raw-aggregate(mean,coarse,fine,Resampling.nearest)
    supported=aggregate(eligible.astype('float32'),coarse,fine,Resampling.nearest)==1
    centered[~supported]=np.nan
    # Units are dimensionless: standardize by scene-wide valid contrast variability.
    sd=np.nanstd(centered)
    return centered/sd if np.isfinite(sd) and sd>0 else np.full_like(centered,np.nan)
