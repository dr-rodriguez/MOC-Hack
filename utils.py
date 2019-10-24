# Utility functions
import pandas as pd
import numpy as np
from astropy.coordinates import ICRS, Galactic, BaseCoordinateFrame
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes.patches import _rotate_polygon
import astropy.units as u
from mocpy import core
import cdshealpix
try:
    from astropy_healpix import HEALPix
except ImportError:
    pass


def query_db(conn, query):
    return pd.read_sql(query, conn)


def convert_to_polygon(center_ra, center_dec, radius, resolution=16):
    """
    Convert a circle to a polygon

    :param center_ra: astropy.units.Quantity
    :param center_dec: astropy.units.Quantity
    :param radius: astropy.units.Quantity
    :param resolution: int
    :return:
    """
    lon = np.linspace(0., 2 * np.pi, resolution + 1)[:-1] * u.radian
    lat = np.repeat(0.5 * np.pi - radius.to_value(u.radian), resolution) * u.radian
    lon, lat = _rotate_polygon(lon, lat, center_ra, center_dec)
    lon = lon.to_value(u.deg).tolist()
    lat = lat.to_value(u.deg).tolist()
    return lon, lat


def parse_s_region(s_region):
    ra = []
    dec = []
    counter = 0

    if s_region is None or s_region.split()[0].upper() not in ('POLYGON','CIRCLE'):
        print('Unsupported shape')
        return None

    if s_region.split()[0].upper() == 'POLYGON':
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
    elif s_region.split()[0].upper() == 'CIRCLE':
        center_ra, center_dec, radius = None, None, None
        for elem in s_region.strip().split():
            try:
                value = float(elem)
            except ValueError:
                continue
            if counter % 2 == 1:
                center_dec = value
            if center_ra is None and counter % 2 == 0:
                center_ra = value
            else:
                radius = value
            counter += 1
        ra, dec = convert_to_polygon(center_ra*u.deg, center_dec*u.deg, radius*u.deg)

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
