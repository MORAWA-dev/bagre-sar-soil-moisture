"""Earth Engine collection setup; imported only in authenticated notebooks."""
import ee
import pandas as pd
import sar_workflow as sw

def initialize(project):
    try: ee.Initialize(project=project)
    except ee.EEException:
        ee.Authenticate(); ee.Initialize(project=project)
    print('Connected:',ee.Number(1).getInfo())

def region(bounds=sw.REGION): return ee.Geometry.Rectangle(list(bounds),proj='EPSG:4326',geodesic=False)
def projection(g): return ee.Projection(g['crs'],g['transform'])

def radar_collection():
    return (ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(region()).filterDate('2024-01-01','2026-01-01')
      .filter(ee.Filter.eq('instrumentMode','IW')).filter(ee.Filter.eq('resolution_meters',10))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH')))

def optical_collection():
    return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(region()).filterDate('2024-01-01','2026-01-01')

def moisture_collection():
    return ee.ImageCollection('NASA/SMAP/SPL3SMP_E/006').filterBounds(region()).filterDate('2024-01-01','2026-01-01')

def track_collection(pass_name,orbit):
    return radar_collection().filter(ee.Filter.eq('orbitProperties_pass',pass_name)).filter(ee.Filter.eq('relativeOrbitNumber_start',int(orbit)))

def track_audit():
    """Choose a pass and relative orbit before examining soil-moisture labels."""
    pilot=region(sw.PILOT); col=radar_collection().filterBounds(pilot)
    def row(im):
        return im.set('audit_row',ee.List([im.get('orbitProperties_pass'),im.get('relativeOrbitNumber_start'),im.date().format('YYYY-MM')]))
    rows=col.map(row).aggregate_array('audit_row').getInfo()
    assert rows, 'No dual-pol IW acquisitions intersect the pilot.'
    df=pd.DataFrame(rows,columns=['pass','orbit','month'])
    out=df.groupby(['pass','orbit']).agg(months=('month','nunique'),scenes=('month','size')).reset_index()
    # Estimate regional coverage on the union of the first month's footprints for each track.
    cov=[]
    aoi=region(); area=aoi.area(100).getInfo()
    for p,o in out[['pass','orbit']].itertuples(index=False,name=None):
        track=track_collection(p,o); first=ee.Image(track.sort('system:time_start').first()).date()
        footprint=track.filterDate(first.advance(-1,'day'),first.advance(14,'day')).geometry(1000)
        cov.append(float(footprint.intersection(aoi,1000).area(1000).getInfo()/area))
    out['footprint_fraction']=cov
    return out.sort_values(['months','footprint_fraction','scenes'],ascending=False).reset_index(drop=True)

def reduce_fine(im,src,dst):
    """Aggregate a single image in its native projection before monthly averaging."""
    im=im.setDefaultProjection(src)
    return im.reduceResolution(reducer=ee.Reducer.mean(),maxPixels=4096).reproject(dst)

def sar_prepare(im,g):
    p=im.select('VV').projection(); dst=projection(g)
    vv=im.select('VV'); vh=im.select('VH'); angle=im.select('angle').reproject(p)
    valid=vv.gt(-35).And(vh.gt(-43)).And(angle.gte(29)).And(angle.lte(46))
    # Convert dB to linear power BEFORE any temporal/spatial mean.
    vvlin=vv.multiply(0.2302585092994046).exp()
    vhlin=vh.multiply(0.2302585092994046).exp()
    stack=vvlin.rename('vv_power').addBands(vhlin.rename('vh_power')).addBands(angle.rename('angle')).updateMask(valid)
    support=reduce_fine(stack.select('vv_power').mask().unmask(0,False),p,dst)
    return reduce_fine(stack,p,dst).updateMask(support.gte(.5)).copyProperties(im,['system:time_start'])

def optical_prepare(im,g):
    p=im.select('B8').projection(); dst=projection(g); scl=im.select('SCL').reproject(p)
    good=scl.eq(4).Or(scl.eq(5)) # vegetation and bare ground only; exclude clouds, water, shadows
    ndvi=im.normalizedDifference(['B8','B4']).rename('ndvi').updateMask(good)
    support=reduce_fine(ndvi.mask().unmask(0,False),p,dst)
    return reduce_fine(ndvi,p,dst).updateMask(support.gte(.5)).copyProperties(im,['system:time_start'])

def smap_prepare(im):
    raw=im.select('soil_moisture_am'); good=im.select('retrieval_qual_flag_am').eq(0)
    return raw.updateMask(good.And(raw.gte(0)).And(raw.lte(.65))).rename('sm').copyProperties(im,['system:time_start'])

def monthly(col,names,count_name,proj):
    empty=ee.Image.constant([0]*len(names)).rename(names).toFloat().updateMask(ee.Image.constant(0)).setDefaultProjection(proj)
    # Normalize every real image to the same unrestricted float band type as
    # the masked fallback. NDVI otherwise has a bounded Float<-1,1> type.
    col=ee.ImageCollection(col).map(lambda im: im.select(names).toFloat()).merge(ee.ImageCollection.fromImages([empty]))
    mean=col.select(names).mean().setDefaultProjection(proj)
    count=col.select(names[0]).count().rename(count_name).unmask(0).setDefaultProjection(proj)
    return mean.addBands(count)

def month_bounds(month):
    start=pd.Timestamp(month); return start.strftime('%Y-%m-%d'),(start+pd.offsets.MonthBegin(1)).strftime('%Y-%m-%d')

def month_images(month,conf):
    start,end=month_bounds(month); fine=conf['fine']; coarse=conf['coarse']
    radar=track_collection(conf['pass'],conf['orbit']).filterDate(start,end).map(lambda im:sar_prepare(im,fine))
    optical=optical_collection().filterDate(start,end).map(lambda im:optical_prepare(im,fine))
    smap=moisture_collection().filterDate(start,end).map(smap_prepare)
    return (monthly(radar,sw.S1[:3],'sar_count',projection(fine)),monthly(optical,['ndvi'],'optical_count',projection(fine)),monthly(smap,['sm'],'smap_days',projection(coarse)))

def land_image(g):
    lc=ee.ImageCollection('ESA/WorldCover/v200').first().select('Map'); p=lc.projection()
    land=lc.neq(80).rename('land_fraction'); crop=lc.eq(40).rename('crop_fraction')
    return reduce_fine(land.addBands(crop),p,projection(g))

def acquire_month(root,month,conf):
    key=str(month)[:7]; a,b,c=month_images(month,conf)
    sw.download(a,root/'data'/f'sar_{key}.tif',conf['fine'],sw.S1)
    sw.download(b,root/'data'/f'optical_{key}.tif',conf['fine'],sw.S2)
    sw.download(c,root/'data'/f'smap_{key}.tif',conf['coarse'],['sm','smap_days'])
