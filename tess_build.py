# Script to build MOC, one per sector
import matplotlib
matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import os
import time
from mocpy import MOC, WCS
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from secrets import CAOM_OPS
from utils import *
from matplotlib import pyplot as plt
from concurrent.futures import ThreadPoolExecutor
import re
plt.interactive(False)

MAX_DEPTH = 9


def get_polygon_moc(row):
    print('{} ({} {}) {}'.format(row['obs_id'], row['s_ra'], row['s_dec'], row['s_region']))
    lon = u.Quantity(row['coords']['ra'] * u.deg)
    lat = u.Quantity(row['coords']['dec'] * u.deg)
    temp_moc = MOC.from_polygon(lon, lat, max_depth=MAX_DEPTH)

    return temp_moc


def get_data(sector):
    query = "SELECT top 1000 * FROM obsPointing WHERE obs_collection='TESS' AND dataproduct_type='image' " \
            "AND sequence_number={}".format(sector)
    df = query_db(CAOM_OPS, query)

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

    return moc


# Build the MOCs
moc_list = []
for s in range(17,19):
    print(s)
    moc = get_data(s)
    moc.write('data/tess_S{}_{}.fits'.format(s, MAX_DEPTH), format='fits', overwrite=True)
    moc_list.append(moc)


# Load MOCs
moc_list = []
for filename in os.listdir('data'):
    if filename.endswith('fits'):
        t = re.findall(r'tess_S(\d+)_', filename)
        s = t[0]
        if int(s) < 17: continue
        moc = MOC.from_fits(os.path.join('data', filename))
        moc_list.append((moc,s))


# Generate the figure(s)
for t in moc_list:
    moc, s = t
    print(s)

    my_plot(moc, frame=Galactic(), save='figures/tess_S{:04d}.png'.format(int(s)))


    # fig = plt.figure(111, figsize=(12.5, 10.5))
    #
    # with WCS(fig,
    #          fov=360 * u.deg,
    #          center=SkyCoord(0, 0, unit='deg', frame='galactic'),
    #          coordsys="galactic",
    #          rotation=Angle(0, u.degree),
    #          projection="AIT") as wcs:
    #     ax = fig.add_subplot(1, 1, 1, projection=wcs)
    #
    #     # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    #     moc.fill(ax=ax, wcs=wcs, alpha=1, linewidth=1.5, fill=True, color="black")
    #     # moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    #
    # plt.tight_layout()
    # plt.savefig('figures/tess_S{:04d}.png'.format(int(s)))


union_moc = MOC.union(*[s[0] for s in moc_list])
my_plot(union_moc, frame=Galactic(), save='caom_TESS_galactic.png')

