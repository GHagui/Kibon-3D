"""Reduz os rasters originais a versões leves (Cloud Optimized GeoTIFF) e completa o shapefile.

Uso: python tools/otimizar_dados.py <pasta_originais> <pasta_saida>
  MDT_kibon.tif (~950 MB, 7 cm)        -> MDT_kibon.tif   (0,5 m, float32, DEFLATE)
  Ortofoto_Kibon*.tif (~4,5 GB, 3 cm)  -> Ortofoto_Kibon.tif (0,2 m, JPEG + máscara)
  Curvas_1m.shp/.dbf (sem .shx/.prj)   -> Curvas_1m.shp/.shx/.dbf/.prj (SIRGAS 2000 / UTM 23S)
"""
import sys, glob, os
import numpy as np, rasterio, shapefile
from rasterio.enums import Resampling
from rasterio.io import MemoryFile
from rasterio.shutil import copy as rio_copy
from pyproj import CRS

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)

def resample(path, res_m, resampling):
    with rasterio.open(path) as d:
        # rasters are in EPSG:4326; convert target metres to degrees at the site latitude
        lat = (d.bounds.top + d.bounds.bottom) / 2
        rx = res_m / (111320 * np.cos(np.radians(lat))); ry = res_m / 110540
        w = round((d.bounds.right - d.bounds.left) / rx); h = round((d.bounds.top - d.bounds.bottom) / ry)
        arr = d.read(out_shape=(d.count, h, w), resampling=resampling)
        prof = d.profile.copy()
        prof.update(width=w, height=h, transform=d.transform * d.transform.scale(d.width / w, d.height / h), driver='GTiff', tiled=True)
        return arr, prof

def write_cog(arr, prof, dst, **opts):
    with MemoryFile() as mem:
        with mem.open(**prof) as tmp:
            tmp.write(arr)
            if prof['count'] == 4:
                from rasterio.enums import ColorInterp as CI
                tmp.colorinterp = [CI.red, CI.green, CI.blue, CI.alpha]
        with mem.open() as tmp:
            rio_copy(tmp, dst, driver='COG', **opts)
    print(dst, os.path.getsize(dst) // 1024 // 1024, 'MB')

arr, prof = resample(os.path.join(src, 'MDT_kibon.tif'), 0.5, Resampling.average)
write_cog(arr, prof, os.path.join(out, 'MDT_kibon.tif'), COMPRESS='DEFLATE', PREDICTOR='YES', LEVEL=9, RESAMPLING='AVERAGE')

orto = sorted(glob.glob(os.path.join(src, 'Ortofoto_Kibon*.tif')))[0]
arr, prof = resample(orto, 0.2, Resampling.average)
write_cog(arr, prof, os.path.join(out, 'Ortofoto_Kibon.tif'), COMPRESS='JPEG', QUALITY=85, RESAMPLING='AVERAGE')

r = shapefile.Reader(shp=open(os.path.join(src, 'Curvas_1m.shp'), 'rb'), dbf=open(os.path.join(src, 'Curvas_1m.dbf'), 'rb'))
with shapefile.Writer(os.path.join(out, 'Curvas_1m'), shapeType=r.shapeType) as w:
    w.fields = r.fields[1:]
    for shp, rec in zip(r.iterShapes(), r.iterRecords()):
        w.shape(shp); w.record(*rec)
with open(os.path.join(out, 'Curvas_1m.prj'), 'w') as f:
    f.write(CRS.from_epsg(31983).to_wkt('WKT1_ESRI'))
print('curvas', len(r))
