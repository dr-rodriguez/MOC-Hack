# Utility functions
import pandas as pd
import numpy as np
from astropy.coordinates import ICRS, Galactic, BaseCoordinateFrame
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import core
import cdshealpix
try:
    from astropy_healpix import HEALPix
except ImportError:
    pass


def query_db(conn, query):
    return pd.read_sql(query, conn)


def parse_s_region_polygon(s_region):
    ra = []
    dec = []
    counter = 0

    if s_region is None or s_region.split()[0].upper() != 'POLYGON':
        print('Not a polygon')
        return None

    for elem in s_region.strip().split():
        try:
            value = float(elem)
        except ValueError:
            continue
        if counter % 2 == 0:
            ra.append(value)
        else:
            dec.append(value)
        counter += 1

    return {'ra': ra, 'dec': dec}


def my_plot(moc, frame=None, labels=False, title='', grid=False, save=''):
    frame = Galactic() if frame is None else frame

    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.pyplot as plt

    plot_order = 8
    if moc.max_order > plot_order:
        plotted_moc = moc.degrade_to_order(plot_order)
    else:
        plotted_moc = moc

    num_pixels_map = 1024
    delta = 2. * np.pi / num_pixels_map

    x = np.arange(-np.pi, np.pi, delta)
    y = np.arange(-np.pi / 2, np.pi / 2, delta)
    lon_rad, lat_rad = np.meshgrid(x, y)
    hp = HEALPix(nside=(1 << plotted_moc.max_order), order='nested')

    if frame and not isinstance(frame, BaseCoordinateFrame):
        raise ValueError("Only Galactic/ICRS coordinate systems are supported."
                         "Please set `coord` to either 'C' or 'G'.")

    pix_map = hp.lonlat_to_healpix(lon_rad * u.rad, lat_rad * u.rad)

    m = np.zeros(12 * 4 ** (plotted_moc.max_order))
    pix_id = core.flatten_pixels(plotted_moc._interval_set._intervals, plotted_moc.max_order)

    # change the HEALPix cells if the frame of the MOC is not the same as the one associated with the plot method.
    if isinstance(frame, Galactic):
        lon, lat = hp.boundaries_lonlat(pix_id, step=2)
        sky_crd = SkyCoord(lon, lat, unit='deg')
        pix_id = hp.lonlat_to_healpix(sky_crd.galactic.l, sky_crd.galactic.b)

    m[pix_id] = 1

    z = np.flip(m[pix_map], axis=1)

    figsize = (12, 10)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="aitoff")

    color_map = LinearSegmentedColormap.from_list('w2r', ['white', 'black'])
    color_map.set_under('w')
    color_map.set_bad('w')

    # Note I flip x
    ax.pcolormesh(x, y, z, cmap=color_map, vmin=0, vmax=1)
    if labels: ax.tick_params(labelsize=14, labelcolor='#000000')
    if title: plt.title(title)
    if grid: plt.grid(True, linestyle='--', linewidth=1, color='#555555')

    ax.set_xticklabels(['210', '240', '270', '300', '330', '0', '30', '60', '90', '120', '150'])
    ax.grid(grid)

    plt.tight_layout()

    if not labels:  # disable tick labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    if save:
        plt.savefig(save)
    else:
        plt.show()
