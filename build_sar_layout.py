"""Portable QGIS seasonal comparison, using actual exported GeoTIFFs."""
import os,sys,shutil
from pathlib import Path
mac=Path('/Applications/QGIS.app/Contents/Resources/qgis')
if mac.is_dir(): os.environ.setdefault('PROJ_DATA',str(mac/'proj'))
from qgis.core import *
from qgis.PyQt.QtGui import QColor,QFont
from qgis.PyQt.QtCore import Qt
source,out=map(lambda x:Path(x).resolve(),sys.argv[1:3]); out.mkdir(parents=True,exist_ok=True)
rasters=out/'rasters'; rasters.mkdir(exist_ok=True)
app=QgsApplication([],False)
if mac.is_dir(): app.setPrefixPath(str(mac),True)
app.initQgis(); project=QgsProject.instance(); crs=QgsCoordinateReferenceSystem('EPSG:32630'); project.setCrs(crs)
wet=['#ffffd9','#c7e9b4','#7fcdbb','#41b6c4','#225ea8','#081d58']; div=['#a50026','#f46d43','#fddbc7','#f7f7f7','#d1e5f0','#4393c3','#053061']
styles={'estimated_sm':(0,.35,wet,'Seasonal estimate (m³/m³)'),'smap_paired':(0,.35,wet,'Paired SMAP (m³/m³)'),'difference':(-.15,.15,div,'Estimate minus SMAP (m³/m³)'),'relative_contrast':(-2,2,div,'Relative score (unitless)'),'coverage':(0,1,wet,'Valid coverage'),'vv_db':(-25,0,wet,'VV (dB)'),'ndvi':(0,.8,wet,'NDVI')}
layers={}
for p in sorted(source.glob('*.tif')):
    key,date=p.stem.rsplit('_',1); lo,hi,colors,title=styles[key]; dst=rasters/p.name; shutil.copy2(p,dst)
    l=QgsRasterLayer(str(dst),f'{date} | {title}'); assert l.isValid(),p
    sh=QgsColorRampShader(lo,hi); sh.setColorRampType(QgsColorRampShader.Interpolated); sh.setColorRampItemList([QgsColorRampShader.ColorRampItem(lo+(hi-lo)*i/(len(colors)-1),QColor(c),'') for i,c in enumerate(colors)])
    rs=QgsRasterShader(lo,hi); rs.setRasterShaderFunction(sh); r=QgsSingleBandPseudoColorRenderer(l.dataProvider(),1,rs); r.setClassificationMin(lo); r.setClassificationMax(hi); l.setRenderer(r)
    project.addMapLayer(l); layers[key,date]=l; project.layerTreeRoot().findLayer(l.id()).setItemVisibilityChecked(False)
extent=layers['relative_contrast','2025-01'].extent(); project.viewSettings().setDefaultViewExtent(QgsReferencedRectangle(extent,crs)); project.layerTreeRoot().findLayer(layers['relative_contrast','2025-01'].id()).setItemVisibilityChecked(True)
layout=QgsPrintLayout(project); layout.initializeDefaults(); layout.setName('Seasonal soil moisture comparison'); layout.pageCollection().pages()[0].setPageSize('A3',QgsLayoutItemPage.Landscape)
def place(it,x,y,w,h):
    layout.addLayoutItem(it);it.attemptMove(QgsLayoutPoint(x,y));it.attemptResize(QgsLayoutSize(w,h));return it
def label(s,x,y,w,h,size=10,bold=False):
    it=QgsLayoutItemLabel(layout);it.setText(s);it.setMargin(0);f=QFont('Helvetica Neue');f.setPointSizeF(size);f.setBold(bold);it.setFont(f);it.setFontColor(QColor('#173c42'));it.setVAlign(Qt.AlignVCenter);return place(it,x,y,w,h)
def box(x,y,w,h,c):
    it=QgsLayoutItemShape(layout);it.setShapeType(QgsLayoutItemShape.Rectangle);it.setSymbol(QgsFillSymbol.createSimple({'color':c,'outline_style':'no'}));return place(it,x,y,w,h)
def ramp(x,y,colors,lo,hi):
    for i in range(60):
        t=i/59*(len(colors)-1);j=min(int(t),len(colors)-2);f=t-j;a,b=QColor(colors[j]),QColor(colors[j+1]);c=QColor(*[round(u+(v-u)*f) for u,v in zip(a.getRgb()[:3],b.getRgb()[:3])]);box(x+i*1.25,y,1.26,3,c.name())
    label(str(lo),x,y+3,20,5,9);label(str(hi),x+62,y+3,22,5,9)
label('Soil moisture and seasonal variability | Bagré, Burkina Faso',10,5,400,12,23,True)
label('January and June 2025 | Sentinel-1 / Sentinel-2 / SMAP | WGS 84 / UTM 30N',10,17,400,8,12)
keys=['estimated_sm','smap_paired','difference','relative_contrast']
for i,key in enumerate(keys):
    x=28+i*96;lo,hi,colors,title=styles[key];label(title,x,30,94,8,11,True);ramp(x,40,colors,lo,hi)
    for j,date in enumerate(['2025-01','2025-06']):
        y=53+j*105;it=QgsLayoutItemMap(layout);place(it,x,y,88,92);it.setCrs(crs);it.setLayers([layers[key,date]]);it.setKeepLayerSet(True);it.storeCurrentLayerStyles();it.setKeepLayerStyles(True);it.zoomToExtent(extent);it.setFrameEnabled(True);it.setFrameStrokeColor(QColor('#aabbbb'))
        if i==0 and j==0: reference=it
label('JAN',7,88,20,10,13,True);label('JUN',7,193,20,10,13,True)
label('↑ N (grid north)',10,255,70,7,10,True);scale=QgsLayoutItemScaleBar(layout);scale.setStyle('Single Box');scale.setLinkedMap(reference);scale.applyDefaultSize();scale.setUnits(Qgis.DistanceUnit.Kilometers);scale.setUnitsPerSegment(20);scale.setNumberOfSegments(2);scale.setNumberOfSegmentsLeft(0);scale.setUnitLabel('km');scale.setHeight(2);layout.addLayoutItem(scale);scale.attemptMove(QgsLayoutPoint(10,264))
label('85 / 79 mapped coarse cells | 70 / 64 paired SMAP cells | 50,032 / 46,169 crop-score pixels (Jan / Jun)',95,253,310,7,11,True)
label('Seasonal baseline is spatially constant within a month. White = masked or unsupported.',95,261,310,7,10)
label('250 m contrast: exploratory, dimensionless, crop fraction ≥50%; not validated soil moisture.',95,269,310,7,10)
label('Ibrahim Kabore | November 2026 | SMAP enhanced 9 km grid; no independent field validation.',10,285,400,7,10)
layout.setReferenceMap(reference);project.layoutManager().addLayout(layout);project.setFilePathStorage(Qgis.FilePathType.Relative);project.setFileName(str(out/'Bagre_SAR_2025.qgz'));assert project.write()
export=QgsLayoutExporter(layout);img=QgsLayoutExporter.ImageExportSettings();img.dpi=300;pdf=QgsLayoutExporter.PdfExportSettings();pdf.dpi=300
assert export.exportToImage(str(out/'Bagre_SAR_2025.png'),img)==QgsLayoutExporter.Success
assert export.exportToPdf(str(out/'Bagre_SAR_2025.pdf'),pdf)==QgsLayoutExporter.Success
print('Exported',len(layers),'layers, QGIS project, PNG and PDF')
