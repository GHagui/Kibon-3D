"""Gera os arquivos do visualizador web (assets/) a partir de dados/.

Uso: python tools/build_assets.py [dados] [assets]

Saídas:
  meta.json        extensão, alturas, índice das curvas, offsets dos binários
  grid.bin.gz      MDT 1024 px (uint16 cm, delta por linha) + direção de escoamento D8 + lâmina de poças
  overlays.bin.gz  curvas de nível, postes (OSM) e luzes das casas
  orto-4k.jpg / orto-2k.jpg  ortofoto recortada no MDT
"""
import sys, os, io, json, gzip, math, heapq
import numpy as np, rasterio, shapefile
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from pyproj import Transformer
from PIL import Image, ImageDraw, ImageFilter

DADOS = sys.argv[1] if len(sys.argv) > 1 else 'dados'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'assets'
GRID = 1024
os.makedirs(OUT, exist_ok=True)

# ---------- terrain ----------
with rasterio.open(f'{DADOS}/MDT_kibon.tif') as d:
    b = d.bounds
    rows = round(GRID * d.height / d.width)
    dem = d.read(1, out_shape=(rows, GRID), resampling=Resampling.average).astype(np.float64)
cols = GRID
valid = (dem > -1000) & np.isfinite(dem)
zmin, zmax = float(dem[valid].min()), float(dem[valid].max())
lon0, lat0 = (b.left + b.right) / 2, (b.top + b.bottom) / 2
R = 6378137.0
kx = math.radians(1) * R * math.cos(math.radians(lat0)); ky = math.radians(1) * R
width, height = (b.right - b.left) * kx, (b.top - b.bottom) * ky
q = np.zeros(dem.shape, np.uint16)
q[valid] = np.clip(np.round((dem[valid] - zmin) * 100) + 1, 1, 65535).astype(np.uint16)

# ---------- hydrology: priority-flood fill, each cell drains to the neighbour it was reached from ----------
NB = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
h = np.where(valid, (q.astype(np.float64) - 1) / 100, np.nan)
fill = np.where(valid, h, 0.0)
dirc = np.full(dem.shape, 255, np.uint8)
seen = ~valid.copy()
pad = np.pad(valid, 1, constant_values=False)
edge = valid & ~(pad[:-2, :-2] & pad[:-2, 1:-1] & pad[:-2, 2:] & pad[1:-1, :-2] & pad[1:-1, 2:] & pad[2:, :-2] & pad[2:, 1:-1] & pad[2:, 2:])
heap = [(h[j, i], j, i) for j, i in zip(*np.nonzero(edge))]
heapq.heapify(heap); seen |= edge
fl, sn, hh, dc = fill, seen, h, dirc
while heap:
    f, j, i = heapq.heappop(heap)
    for n, (a, c) in enumerate(NB):
        nj, ni = j + c, i + a
        if 0 <= nj < rows and 0 <= ni < cols and not sn[nj, ni]:
            sn[nj, ni] = True
            dc[nj, ni] = 7 - n          # neighbour drains back towards (j, i)
            v = max(hh[nj, ni], f + 1e-4)
            fl[nj, ni] = v
            heapq.heappush(heap, (v, nj, ni))
pond = np.zeros(dem.shape, np.uint8)
pond[valid] = np.clip(np.round((fill[valid] - h[valid]) * 100), 0, 255).astype(np.uint8)
delta = np.diff(q.astype(np.int32), axis=1, prepend=0).astype(np.int16)   # compresses far better than raw heights
grid_bin = delta.tobytes() + dirc.tobytes() + pond.tobytes()
open(f'{OUT}/grid.bin.gz', 'wb').write(gzip.compress(grid_bin, 9))

# ---------- orthophoto cropped to the DEM ----------
with rasterio.open(f'{DADOS}/Ortofoto_Kibon.tif') as o:
    win = from_bounds(b.left, b.bottom, b.right, b.top, o.transform)
    th = round(4096 * rows / cols)
    rgb = o.read([1, 2, 3], out_shape=(3, th, 4096), window=win, resampling=Resampling.average, boundless=True, fill_value=0)
    msk = o.read_masks(1, out_shape=(th, 4096), window=win, resampling=Resampling.average, boundless=True)
img = np.transpose(rgb, (1, 2, 0)).copy()
img[msk < 128] = (24, 26, 22)
ortho = Image.fromarray(img)
for size, name in ((4096, 'orto-4k.jpg'), (2048, 'orto-2k.jpg')):
    # square power-of-two textures stay under the 4096 limit of most phones
    ortho.resize((size, size), Image.LANCZOS).save(f'{OUT}/{name}', 'JPEG', quality=84, optimize=True, progressive=True)

# ---------- contours (shapefile in SIRGAS 2000 / UTM 23S) ----------
tr = Transformer.from_crs(31983, 4326, always_xy=True)
to_local = lambda lon, lat: ((np.asarray(lon) - lon0) * kx, -(np.asarray(lat) - lat0) * ky)
r = shapefile.Reader(f'{DADOS}/Curvas_1m')
cpts, lines = [], []
for shp, rec in zip(r.iterShapes(), r.iterRecords()):
    for n, start in enumerate(shp.parts):
        end = shp.parts[n + 1] if n + 1 < len(shp.parts) else len(shp.points)
        p = np.array(shp.points[start:end])
        if len(p) < 2: continue
        x, z = to_local(*tr.transform(p[:, 0], p[:, 1]))
        lines.append([sum(len(c) for c in cpts) // 2, len(p), float(rec['ELEVATION'])])
        cpts.append(np.column_stack([x, z]).ravel())
contours = np.round(np.concatenate(cpts) * 20).astype(np.int16)        # 5 cm, x east / z south

# ---------- street lights along OSM roads ----------
osm = json.load(open(f'{DADOS}/osm_vias.json'))
KIND = {'primary': (24, 2), 'primary_link': (24, 2), 'secondary': (25, 2), 'tertiary': (27, 2),
        'residential': (30, 1), 'unclassified': (30, 1), 'living_street': (30, 1), 'service': (38, 1),
        'footway': (34, 0), 'steps': (20, 0), 'pedestrian': (30, 0), 'path': (45, 0)}
def survey(x, z):
    i, j = round((x / width + 0.5) * (cols - 1)), round((z / height + 0.5) * (rows - 1))
    return 0 <= i < cols and 0 <= j < rows and valid[j, i]
lx, lk = [], []
for e in osm['elements']:
    kind = KIND.get(e.get('tags', {}).get('highway'))
    if not kind or 'geometry' not in e: continue
    step, cls = kind
    pts = [tuple(map(float, to_local(g['lon'], g['lat']))) for g in e['geometry']]
    carry, side = step / 2, 1
    for (x0, z0), (x1, z1) in zip(pts, pts[1:]):
        L = math.hypot(x1 - x0, z1 - z0)
        if L == 0: continue
        ux, uz = (x1 - x0) / L, (z1 - z0) / L
        s = carry
        while s < L:
            off = 3.5 * side if cls else 0.0
            x, z = x0 + ux * s - uz * off, z0 + uz * s + ux * off
            if survey(x, z): lx += [x, z]; lk.append(cls)
            side = -side; s += step
        carry = s - L
lamps = np.round(np.array(lx) * 20).astype(np.int16)

# ---------- window lights: OSM has almost no buildings here, so roofs come from the orthophoto ----------
RES = 0.7
TW, TH = round(width / RES), round(height / RES)
a = np.asarray(ortho.resize((TW, TH), Image.BILINEAR)).astype(np.int16)
rr, gg, bb = a[..., 0], a[..., 1], a[..., 2]
ok = ~((abs(rr - 24) < 6) & (abs(gg - 26) < 6) & (abs(bb - 22) < 6))
veg = (2 * gg - rr - bb > 14) | ((gg > rr) & (gg > bb + 12))
lum = (rr + gg + bb) / 3
built = ok & ~veg & (lum > 45) & (lum < 240)
# roofs are busy (edges, shadows between houses); bare soil, sports fields and water are smooth
blur = lambda arr, rad: np.asarray(Image.fromarray(np.clip(arr, 0, 255).astype(np.uint8)).filter(ImageFilter.BoxBlur(rad)), np.float32)
m1 = blur(lum, 3); m2 = blur(lum ** 2 / 255, 3) * 255
std = blur(np.sqrt(np.maximum(0, m2 - m1 ** 2)) * 4, 5) / 4
built &= std > 14
road = Image.new('L', (TW, TH), 0); dr = ImageDraw.Draw(road)
WID = {'primary': 12, 'primary_link': 10, 'secondary': 11, 'tertiary': 9, 'residential': 7,
       'unclassified': 7, 'service': 5, 'footway': 2.5, 'steps': 2.5, 'path': 2}
for e in osm['elements']:
    w = WID.get(e.get('tags', {}).get('highway'))
    if w and 'geometry' in e:
        px = [((x / width + 0.5) * TW, (z / height + 0.5) * TH) for x, z in (tuple(map(float, to_local(g['lon'], g['lat']))) for g in e['geometry'])]
        dr.line(px, fill=255, width=max(1, round(w / RES)))
built &= np.asarray(road) == 0
dense = built & (blur(built * 255, 6) / 255 > 0.55)
rng = np.random.default_rng(7)
hx, hk = [], []
SP = 6.0
for gz in np.arange(-height / 2, height / 2, SP):
    for gx in np.arange(-width / 2, width / 2, SP):
        x, z = gx + rng.uniform(0, SP), gz + rng.uniform(0, SP)
        px, py = int((x / width + 0.5) * TW), int((z / height + 0.5) * TH)
        if 0 <= px < TW and 0 <= py < TH and dense[py, px] and rng.random() < 0.7 and survey(x, z):
            hx += [x, z]
            hk.append(rng.choice(4, p=[0.55, 0.3, 0.1, 0.05]) * 4 + rng.choice(3, p=[0.5, 0.35, 0.15]))
houses = np.round(np.array(hx) * 20).astype(np.int16)

parts = [contours.tobytes(), lamps.tobytes(), np.array(lk, np.uint8).tobytes(), houses.tobytes(), np.array(hk, np.uint8).tobytes()]
offs, o = [], 0
for p in parts: offs.append(o); o += len(p)
open(f'{OUT}/overlays.bin.gz', 'wb').write(gzip.compress(b''.join(parts), 9))

meta = dict(cols=cols, rows=rows, width=width, height=height, zmin=zmin, zmax=zmax, lon0=lon0, lat0=lat0, kx=kx, ky=ky,
            lines=lines, nContour=len(contours) // 2, nLamps=len(lk), nHouses=len(hk),
            overlays=dict(contours=offs[0], lamps=offs[1], lampKind=offs[2], houses=offs[3], houseKind=offs[4]))
json.dump(meta, open(f'{OUT}/meta.json', 'w'), separators=(',', ':'))
for f in sorted(os.listdir(OUT)): print(f, os.path.getsize(f'{OUT}/{f}') // 1024, 'KB')
print('lamps', len(lk), 'houses', len(hk), 'lines', len(lines))
