# Read the TIC
import time
from mocpy import MOC, WCS
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from secrets import TIC_CONN
from utils import *
from matplotlib import pyplot as plt
from concurrent.futures import ThreadPoolExecutor

MAX_DEPTH = 15


def get_catalog_moc(row):
    lon = u.Quantity(row['ra'] * u.deg)
    lat = u.Quantity(row['dec'] * u.deg)
    temp_moc = MOC.from_lonlat(lon, lat, max_norder=MAX_DEPTH)

    return temp_moc


co = 100000
st = 1
en = st + co
query = "SELECT objID,ra,dec FROM objectCoords WITH (NOLOCK) WHERE objID between {} AND {}".format(st,en)
df = query_db(TIC_CONN, query)

# Generate MOC
start_time = time.time()
pool = ThreadPoolExecutor(max_workers=4)
results = list(pool.map(get_catalog_moc, [row for _, row in df.iterrows()]))
end_time = time.time()
print('Total time : {} seconds'.format(end_time - start_time))

start_time = time.time()
# moc = MOC.union(results[0], *results[1:])
moc = MOC.union(*results)
end_time = time.time()
print('Total time : {} seconds'.format(end_time - start_time))

# Save MOC
# hdulist = moc.serialize(format='fits')
# hdulist.writeto('TIC_v70_{}.fits'.format(MAX_DEPTH))
moc.write('TIC_v70_{}.fits'.format(MAX_DEPTH), format='fits', overwrite=True)

# Read from file
moc = MOC.from_fits('TIC_v70_15.fits')
