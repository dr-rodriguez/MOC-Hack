# Script to build MOC, one per sector
import matplotlib
matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import os
from mocpy import MOC, WCS
import astropy.units as u
from secrets import CAOM_CONN as CAOM
from utils import parse_s_region, my_plot, query_db
from matplotlib import pyplot as plt
from astropy.coordinates import Galactic
plt.interactive(False)

MAX_DEPTH = 9


def get_polygon_moc(row):
    print('{} ({} {}) {}'.format(row['obs_id'], row['s_ra'], row['s_dec'], row['s_region']))
    lon = u.Quantity(row['coords']['ra'] * u.deg)
    lat = u.Quantity(row['coords']['dec'] * u.deg)
    temp_moc = MOC.from_polygon(lon, lat, max_depth=MAX_DEPTH)

    return temp_moc


def get_data():
    query = "SELECT top 1000 * FROM obsPointing WHERE obs_collection='K2' AND dataproduct_type='image'"
    df = query_db(CAOM, query)
    if len(df) == 0:
        return None

    df['coords'] = df.apply(lambda x: parse_s_region(x['s_region']), axis=1)

    # Generate MOC
    results = [get_polygon_moc(row) for _, row in df.iterrows()]

    # Union of MOCs
    if len(results) > 1:
        moc = MOC.union(*results)
    else:
        moc = results

    return moc


# Build the MOCs
moc = get_data()
moc.write('data/k2/k2_{}.fits'.format(MAX_DEPTH), format='fits', overwrite=True)

# Load MOCs
moc = MOC.from_fits(os.path.join('data/k2', 'k2_{}.fits'.format(MAX_DEPTH)))

# Plot MOCs
my_plot(moc, frame=Galactic(), save='caom_K2_galactic.png', grid=True, labels=True, color='blue')

