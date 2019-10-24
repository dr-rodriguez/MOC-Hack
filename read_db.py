import time
from mocpy import MOC, WCS
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
# from secrets import CAOM_CONN
from utils import *
from matplotlib import pyplot as plt
from concurrent.futures import ThreadPoolExecutor

# MAX_DEPTH = 14
MAX_DEPTH = 9
# TESS pixels are ~21" so no need to go lower than level 14 (~13")
# Depth 9 has angular size ~7' This is comparable to HTM level 10


def get_polygon_moc(row):
    print('{} ({} {}) {}'.format(row['obs_id'], row['s_ra'], row['s_dec'], row['s_region']))
    lon = u.Quantity(row['coords']['ra'] * u.deg)
    lat = u.Quantity(row['coords']['dec'] * u.deg)
    temp_moc = MOC.from_polygon(lon, lat, max_depth=MAX_DEPTH) #,
                                # inside=SkyCoord(ra=row['s_ra'] * u.deg, dec=row['s_dec'] * u.deg))

    print(temp_moc.serialize('json'))

    return temp_moc


query = "SELECT top 1000 * FROM obsPointing WHERE obs_collection='TESS' AND dataproduct_type='image'"
# query = "SELECT obs_id, s_ra, s_dec, s_region FROM obsPointing WITH (NOLOCK) WHERE obs_collection IN ('HST','HLA')"
df = query_db(CAOM_CONN, query)

df['coords'] = df.apply(lambda x: parse_s_region(x['s_region']), axis=1)

# Generate MOC
start_time = time.time()
pool = ThreadPoolExecutor(max_workers=4)
results = list(pool.map(get_polygon_moc, [row for _, row in df.iterrows()]))
end_time = time.time()
print('Total time : {} seconds'.format(end_time - start_time))

# Union of MOCs
start_time = time.time()
moc = MOC.union(*results)
end_time = time.time()
print('Total time : {} seconds'.format(end_time - start_time))

# Save MOC
# hdulist = moc.serialize(format='fits')
# hdulist.writeto('tess_s0001_{}.fits'.format(MAX_DEPTH))
moc.write('tess_s0001_{}.fits'.format(MAX_DEPTH), format='fits', overwrite=True)
# moc.write('tess_{}.fits'.format(MAX_DEPTH), format='fits', overwrite=True)
# moc.write('hst_{}.fits'.format(MAX_DEPTH), format='fits', overwrite=True)

# Read from file
moc = MOC.from_fits('tess_s0001_14.fits')
moc = MOC.from_fits('tess_9.fits')


moc.plot(title='TESS')
plt.savefig('full_tess.png')

# Plot the MOC using matplotlib
fig = plt.figure(111, figsize=(10, 8))
# Define a astropy WCS easily
with WCS(fig,
        fov=360 * u.deg,
        center=SkyCoord(0, 0, unit='deg', frame='galactic'),
        coordsys="galactic",
        rotation=Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    moc.fill(ax=ax, wcs=wcs, alpha=1, linewidth=1.5, fill=True, color="green")
    moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")

major_ticks = np.arange(0, 101, 20)
ax.set_xticks(major_ticks)
ax.set_yticks(major_ticks)

# plt.xlabel('RA')
# plt.ylabel('Dec')
# plt.title('Coverage of TESS Sector 1')
plt.grid(color="black", linestyle="dotted")
plt.savefig('temp_full.png')
