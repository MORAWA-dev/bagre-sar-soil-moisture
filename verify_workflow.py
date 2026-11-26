"""Run offline checks with synthetic raster data; no live Earth Engine claims."""
import json,sys, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
ROOT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'code'))
import sar_workflow as sw

for p in sorted((ROOT/'notebooks').glob('*.ipynb')):
    nb=json.loads(p.read_text()); assert nb['nbformat']==4
    for cell in nb['cells']:
        if cell['cell_type']=='code':
            compile('\n'.join(l for l in ''.join(cell['source']).splitlines() if not l.startswith('%')),str(p),'exec')
    print('Compiled',p.name,flush=True)
for p in (ROOT/'code').glob('*.py'): compile(p.read_text(),str(p),'exec')
assert np.allclose(sw.db(sw.power([-20,-10])),[-20,-10])
fine=sw.grid(sw.REGION,sw.CRS,250)
coarse=sw.native_grid({'crs':'EPSG:4326','transform':[.081,0,-180,0,-.081,90]})
print('Synthetic fine and geographic coarse grids:',sw.shape(fine),sw.shape(coarse))
rng=np.random.default_rng(42); sh=sw.shape(fine)
land=np.stack([np.ones(sh),np.full(sh,.7)])
tables=[]; last=None
for m in sw.MONTHS:
    phase=np.sin(2*np.pi*(m.month-1)/12)
    vv=np.clip(.08+.025*phase+rng.normal(0,.02,sh),.005,.3)
    vh=np.clip(.02+.005*phase+rng.normal(0,.006,sh),.001,.1)
    angle=rng.uniform(32,40,sh); ndvi=np.clip(.35+.15*phase+rng.normal(0,.05,sh),0,1)
    sar=np.stack([vv,vh,angle,np.full(sh,3)]); optical=np.stack([ndvi,np.full(sh,4)])
    sm=np.clip(.08+.2*sw.aggregate(ndvi,fine,coarse)+.6*sw.aggregate(vv,fine,coarse)+rng.normal(0,.015,sw.shape(coarse)),0,.65)
    smap=np.stack([sm,np.full(sw.shape(coarse),12)])
    table,valid=sw.summarize_month(sar,optical,land,smap,fine,coarse,m.strftime('%Y-%m')); tables.append(table); last=(sar,optical,valid,table)
df=sw.assign_splits(pd.concat(tables,ignore_index=True)); assert df.eligible.sum()>1000
assert df[df.split=='train'].month.str.startswith('2024').all()
assert df[df.split.isin(['joint_test','temporal_test'])].month.str.startswith('2025').all()
train_cells=set(zip(df[df.split=='train'].row,df[df.split=='train'].col)); test_cells=set(zip(df[df.split=='joint_test'].row,df[df.split=='joint_test'].col)); assert not train_cells&test_cells
assert not df.duplicated(['month','row','col']).any()
print('Split counts:',df[df.eligible].groupby('split').size().to_dict())
models,scores,pred,sel=sw.fit_models(df)
assert len(scores)==24 and scores.n.gt(0).all() and np.isfinite(scores.rmse).all()
print('Six baselines/models × four splits; selected:',sel)
sar,optical,valid,t=last
contrast=sw.fine_contrast(models['fusion_ridge'][1],sar,optical,valid,fine,coarse,t.predictable.to_numpy().reshape(sw.shape(coarse)))
assert contrast.shape==sw.shape(fine) and np.isfinite(contrast).any() and abs(np.nanstd(contrast)-1)<1e-5
with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)/'test.tif'; sw.write_tif(p,contrast,fine,['contrast']); a=sw.read_tif(p,fine,['contrast'])[0]; assert np.allclose(a,contrast,equal_nan=True)
    p=Path(tmp)/'model.joblib'; joblib.dump(models,p); loaded=joblib.load(p); x=df[df.eligible].iloc[:5][sw.FEATURES]; assert np.allclose(loaded['fusion_ridge'][1].predict(x),models['fusion_ridge'][1].predict(x))
print('Raster aggregation, masks, splits, models, persistence and fine contrast passed.')
