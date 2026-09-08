"""

ESO Download Portal
https://www.eso.org/qi/
https://archive.eso.org/downloadportal/0b9bc5e8-9277-4398-bcaf-beba70063390

https://kids.strw.leidenuniv.nl/overview.php

https://zenodo.org/records/14226117

https://astroquery.readthedocs.io/en/latest/api/astroquery.xmatch.XMatchClass.html

https://cds-astro.github.io/tutorials/1_Intro_to_CDS_services_in_notebooks.html

https://tapvizier.cds.unistra.fr/adql/
https://tapvizier.cds.unistra.fr/adql/help.html
https://cds.unistra.fr/help/documentation/vizier-more/adql-vizier/

https://docs.g-vo.org/adql/

"""
# Standard Library
from pathlib import Path
import os
import sys
import time

import argparse
import configparser
import logging
import getpass
import inspect
from inspect import currentframe, getframeinfo
import traceback

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import numpy as np

# Astropy
from astropy import units as u
from astropy.table import Table

# Access astronomical databases
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astroquery.utils.tap.core import TapPlus

import pyvo

import fitsio
from fitsio import FITS,FITSHDR


DEBUG = False
if DEBUG:
    help(XMatch)

""" Example code

get_available_tables(self, *, cache=True)

query(self, cat1, cat2, max_distance, *,
colRA1=None, colDec1=None, colRA2=None, colDec2=None,
area='allsky', cache=True, get_query_payload=False, **kwargs)

query_async(self, cat1, cat2, max_distance, *,
colRA1=None, colDec1=None, colRA2=None, colDec2=None, area='allsky',
cache=True, get_query_payload=False, **kwargs)

"""



def get_logger(name=__name__, level=logging.INFO, test=False):
    """Return a module-level logger, adding a StreamHandler if none exists.

    Parameters
    ----------
    name : str, optional
        Logger name. Defaults to the module name.
    level : int, optional
        Logging level. Defaults to logging.INFO.
    test : bool, optional
        If True, run self-test assertions and demo messages.
        Defaults to False.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(
            fmt='%(asctime)s.%(msecs)03d %(name)-12s %(module)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
            datefmt='%y-%m-%dT%H:%M:%S'
        ))
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.info(f'Logger created: {name}')

    if test:
        logger.debug(f'DEBUG message (should not appear at INFO level)')
        logger.info(f'INFO message')
        logger.warning(f'WARNING message')

        logger2 = logging.getLogger(name)
        assert logger is logger2, 'get_logger is not idempotent'
        assert len(logger.handlers) == 1, (
            f'Expected 1 handler, found {len(logger.handlers)}'
        )
        assert logger.level == level, (
            f'Expected level {level}, got {logger.level}'
        )

        logger.info(f'All tests passed')

        if sys.stdin.isatty():
            input('Enter any key to continue... ')

    return logger


def try_logger(test=False):
    """ check if a logger exists

    """
    import logging

    try:
        logger
    except NameError:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            # logger.addHandler(logging.StreamHandler())
            logger.setLevel(logging.INFO)
        logger.info('Create logger')

    if test:
        logger.info('')

    return logger


def try_logger_v2():
    """Check if a logger exists."""
    import logging

    try:
        logger
    except NameError:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            logger.addHandler(logging.StreamHandler())
            logger.setLevel(logging.INFO)
            logger.propagate = False  # Prevent duplicate messages
        logger.info('Create logger')

    return logger


def kids_make_tile_corners(radec_center=[180.0, 0.0],
                      width=1.0,
                      height=1.0,
                      edgewidth=None):
    """ Warning: probably will not work at poles and if a tile crosses over RA 24->0

    The final calibrated, coadded images from the Astro-WISE pipeline have a
    uniform pixel scale of 0.2 arcsec.

    The r-band detection images, which result from a reduction using THELI,
    specifically designed for optimal image shape measurement and small
    astrometric reductions, are the basis of the multi-band catalogue.
    They have a pixel scale that is closer to the native OmegaCAM pixel
    size of 0.213”.

    An example image has 21000 x 21000 pixels which at 0.2" per pixel is
    4200" x 4200" which is 3600 + 600 or 300" (5') window edge
    70' x 70'

    0.2" per pixel = 5.555555e-5 degrees per pixel

    see:

    https://docs.astropy.org/en/stable/coordinates/matchsep.html

    https://docs.astropy.org/en/stable/wcs/supported_projections.html

    https://docs.astropy.org/en/latest/coordinates/angles.html#wrapping-and-bounds
    https://docs.astropy.org/en/latest/coordinates/angles.html#longitude-and-latitude-objects

    https://docs.astropy.org/en/latest/coordinates/angles.html#with-random-sampling

    """

    import math
    import numpy

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import Angle

    ncorners = 4
    ra_corner = np.zeros(ncorners)
    dec_corner = np.zeros(ncorners)

    center = \
        SkyCoord(radec_center[0]*u.deg, radec_center[1]*u.deg,
                 frame='icrs')

    dec_top = radec_center[1] + (height/2.0)
    dec_bottom = radec_center[1] - (height/2.0)

    # Corner order:
    # Top Left (TL)
    # Top Right (TR)
    # Bottom Right BR)
    # Bottom Left (BL)

    cosdec_correction = (1.0/math.cos(math.radians(dec_top)))
    print(dec_top, cosdec_correction)

    ra_corner[0] = radec_center[0] - (cosdec_correction*width/2.0)
    dec_corner[0] = dec_top

    ra_corner[1] = radec_center[0] + (cosdec_correction*width/2.0)
    dec_corner[1] = dec_top

    cosdec_correction = (1.0/math.cos(math.radians(dec_bottom)))
    print(dec_bottom, cosdec_correction)

    # note order is + -> - to ensure clockwise
    ra_corner[2] = radec_center[0] + (cosdec_correction*width/2.0)
    dec_corner[2] = dec_bottom

    ra_corner[3] = radec_center[0] - (cosdec_correction*width/2.0)
    dec_corner[3] = dec_bottom


    return ra_corner, dec_corner

ra_corner, dec_corner = kids_make_tile_corners(
    radec_center=[180.0, 0.0],
    width=1.0,
    height=1.0,
    edgewidth=None)

# Close the polygon by adding the first point at the end
ra_plot = list(ra_corner) + [ra_corner[0]]
dec_plot = list(dec_corner) + [dec_corner[0]]

fig, ax = plt.subplots()

filled_polygon = True
ax.plot(ra_plot, dec_plot, 'r-', linewidth=2)
if filled_polygon:
    ax.fill(ra_plot, dec_plot, alpha=0.2)



ra_corner, dec_corner = kids_make_tile_corners(
    radec_center=[180.0, 60.0],
    width=1.0,
    height=1.0,
    edgewidth=None)

for i, (ra, dec) in enumerate(zip(ra_corner, dec_corner)):
    print(f"Corner {i}: RA={ra}, Dec={dec}")

 # Close the polygon by adding the first point at the end
ra_plot = list(ra_corner) + [ra_corner[0]]
dec_plot = list(dec_corner) + [dec_corner[0]]

filled_polygon = True
ax.plot(ra_plot, dec_plot, 'r-', linewidth=2)
if filled_polygon:
    ax.fill(ra_plot, dec_plot, alpha=0.2)

ax.set_aspect('equal')

plt.show()

sys.exit()


def tap_get_radec_limits(tap_service_name=None,
                         tap_table=None,
                         colnames_radec=['RAJ2000', 'DEJ2000']):
    """Get the TAP RADec limits using brute force

    https://dp0-2.lsst.io/data-access-analysis-tools/adql-recipes.html

    See: https://kids.strw.leidenuniv.nl/overview.php

    this is quite clumsy but just needs to be run once and the results saved
    Each TAP query took three minutes!

    see: https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html

    https://zenodo.org/records/14226117

    DR5 [all]
    RA range: 0.0 to 360.0
    DEC range: -35.609927 3.961583


    N:
    RA range:
    Dec range: -35.609927 -25.719758

    S:
    RA range:
    Dec range: -35.609927 -25.719758

    """

    t0 = time.time()

    tap_service_name = "http://tapvizier.cds.unistra.fr/TAPVizieR/tap"
    # Connect to VizieR TAP service
    print(f'pyvo.dal.TAPService')
    tap = pyvo.dal.TAPService(tap_service_name)

    # All
    # note \" since CDS TAP needs "tablename"
    query0 = f"""
    SELECT
    -- KIDS All
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (DEJ2000 >= -40.0) AND (DEJ2000 <= 5.0)
    """
    # Result
    #   RA Range: 0.0 360.0
    #   Dec range: -35.609927 3.961583


    # KIDS-N All
    query1 = f"""
    SELECT
    -- KIDS-N All
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (DEJ2000 >= -10.0) AND (DEJ2000 <= 5.0)
    """
    # Result:
    # RA Range: 128.489853 238.510189
    # Dec range: -3.961563 3.961583


    # WAVES region
    # KiDS-N-W2 128.5 to 141.5 -2.0 to +3.0
    query2 = f"""
    SELECT
    -- KiDS-N-W2 WAVES
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (RAJ2000 >= 120.0) AND (RAJ2000 <= 143.0) AND
          (DEJ2000 >= -10.0) AND (DEJ2000 <= 5.0)
    """
    # Result:
    # RA range: 128.489853 141.903234
    # Dec range: -1.983534 2.97255


    # KIDS-N-D2: COSMOS region
    # KiDS-N-D2	149.6 to 151.6	1.7 to 2.7
    query3 = f"""
    SELECT
    -- KIDS-N-D2: COSMOS region
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (RAJ2000 >= 149.0) AND (RAJ2000 <= 152.0) AND
          (DEJ2000 >= -10.0) AND (DEJ2000 <= 5.0)
    """
    # Result:
    # RA Range: 149.572644 150.593402
    # Dec Range: 1.71137 2.700334


    # KIDS-N-D2: COSMOS region
    # KiDS-N-D2	149.6 to 151.6	1.7 to 2.7
    query4 = f"""
    SELECT
    -- KIDS-N:
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (RAJ2000 >= 152.0) AND (RAJ2000 <= 240.0) AND
          (DEJ2000 >= -10.0) AND (DEJ2000 <= 5.0)
    """
    # Result:
    # RA Range:
    # Dec Range:


    # Southern region
    query5 = f"""
    SELECT
    -- KIDS-S All
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (DEJ2000 <= -20.0) AND (DEJ2000 >= -40.0)
    """
    # Result:
    # Ra range: 0.0 360
    # Dec range: -35.609927 -25.719758


    #
    query6 = f"""
    SELECT
    -- KIDS-S > 12hrs
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (DEJ2000 <= -15.0) AND (RAJ2000 >= 180.0)
    """
    # Result: S RA range 328.809452 to 360.0000

    query7 = f"""
    SELECT
    -- KIDS-S < 12hrs
        MIN(RAJ2000), MAX(RAJ2000),
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    WHERE (DEJ2000 <= -15.0) AND (RAJ2000 <= 180.0)
    """
    # Result: S RA range 0 to 54.136983

    # Superset:
    # N:
    # RA range: 128.0 to 239.0
    # Dec raage:

    query8 = f"""
    -- RA range timing test
    SELECT
        MIN(RAJ2000), MAX(RAJ2000)
    FROM \"{tap_table}\"
    """

    # indexing timing test
    query9 = f"""
    SELECT
    -- Dec range timing test
        MIN(DEJ2000), MAX(DEJ2000)
    FROM \"{tap_table}\"
    """
    # Result:
    # N:
    # RA range: 128.0 to 239.0
    # Dec rnage:


    querylist = [query0, query1, query2, query3,
                 query4, query5, query6, query7,
                 query8, query9]

    istart = 3
    iend = 7
    nquerys = len(querylist)
    for iquery, query in enumerate(querylist[istart:iend+1]):

        print(f'query{iquery+istart}: {query}')
        result = tap.search(query)
        print(f'result: {result}')

        print('Elapsed time:', time.time() - t0, 'seconds')

    return


def explore_kids_tile(tap_service_name=None,
                      table_name="II/383/kids_dr5",
                      tile_name='KIDS_150.1_2.2',
                      ntop=None,
                      tap_async_mode=False,
                      test=False):
    """
    explore a single catalogue tile

    """


    return


def tap_get_tile_list(tap_service_name=None,
                      table_name="II/383/kids_dr5",
                      ntop=None,
                      colname='tile',
                      tap_async_mode=False,
                      test=False):
    """
    Note: async is a python reserved word since python 3.7
    https://cdsarc.cds.unistra.fr/viz-bin/cat/II/383

    tile example
    KIDS_1.2_-35.1

    """

    import sqlparse

    test = True
    if test:
        #ntop = 40
        ntop = None
        tap_async_mode = True
        tap_async_mode = False

    logger.info('')
    logger.info(f'TAP service name: {tap_service_name}')
    logger.info(f'TAP async: {tap_async_mode}')
    print(f'Table: {table_name}')


    if tap_service_name is None:
        tap_service_name = "http://tapvizier.cds.unistra.fr/TAPVizieR/tap"
        logger.info(f'TAP service name: {tap_service_name}')

    service = pyvo.dal.TAPService(tap_service_name)

    service.describe()
    print(service.describe)

    if table_name is None:
        table_name = "II/383/kids_dr5"

    # this fails since DISTINCT not supported by TAPVizier
    query = f"""
        SELECT DISTINCT {colname}
            FROM \"{table_name};\"
    """

    if ntop is None:
        query = \
            f"""SELECT tile, COUNT(tile) AS NSources
                FROM \"{table_name}\"
                GROUP BY tile;"""

    if ntop is not None:
        query = \
            f"""SELECT TOP {ntop}  tile, COUNT(tile) AS NSources
                FROM \"{table_name}\"
                GROUP BY tile;"""

    radec_limits = True
    print(f'ntop: {ntop}')
    print(f'radec_limits: {radec_limits}')
    if radec_limits:
        # build query as a list of strings
        query = ["SELECT "] # note the trailing space

        if ntop is not None:
            query.append(f"TOP {ntop} ")


        #query.append(f"""
        #    tile, COUNT(tile) AS NSources,
        #    MIN(RAJ2000) AS RA_min, MAX(RAJ2000) AS RA_max,
        #    MIN(DEJ2000) AS Dec_min, MAX(DEJ2000) AS Dec_max,
        #    MIN(Xpos) AS Xpos_min, MAX(Xpos) AS Xpos_max,
        #    MIN(Ypos) AS Ypos_min, MAX(Ypos) AS Ypos_max
        #    FROM \"{table_name}\"
        #    GROUP BY tile;""")

        query.append(f"""
            tile, COUNT(tile) AS NSources,
            MIN(RAJ2000) AS RA_min, MAX(RAJ2000) AS RA_max,
            MIN(DEJ2000) AS Dec_min, MAX(DEJ2000) AS Dec_max
            FROM \"{table_name}\" """)

        dec_limit = None
        if dec_limit is not None:
            query.append(f"WHERE DEJ2000 > -10.0 ")

        query.append(f"GROUP BY tile;")


        # create a single string query with line breaks from the list of strings
        query = "\n".join(query)


    """
    ORDER BY tile;
    """

    print(f'Tap query: {query}')
    print()

    # create a human readable version of the sql command
    pretty_query = \
        sqlparse.format(query, reindent=True, keyword_case='upper')
    print(pretty_query)
    print()

    query = pretty_query
    t0 = time.time()

    logger.info(f'tap_async_mode: {tap_async_mode}')
    logger.info(f'execute the TAP query')
    if not tap_async_mode:
        result = service.search(query)

    if tap_async_mode:
        job = service.submit_job(query)
        print(job.url)          # job URL on the server
        print(job.phase)

        job.run()
        while job.phase not in ("COMPLETED", "ERROR", "ABORTED"):
            elapsed = time.time() - t0
            print(f'{job.phase} Elapsed time: {elapsed:.2f} seconds')
            time.sleep(10)
            job.raise_if_error()

        job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=3600)
        job.raise_if_error()
        result = job.fetch_result()


    print('Elapsed time:', time.time() - t0, 'seconds')
    num_rows = len(result)
    print(f"Number of rows: {num_rows}")
    nsources = result['NSources']
    print(min(nsources), max(nsources))
    print()

    print('result: ', result)

    nrows = len(result)
    ra = np.zeros(nrows)
    dec = np.zeros(nrows)
    #nsources = np.zeros(nrows)

    for i, row in enumerate(result):
        radec = row[colname].split("_")
        ra[i] = float(radec[1])
        dec[i] = float(radec[2])
        print(f'ra, dec: {i}, {ra[i]}, {dec[i]}')

    tile_list = np.array(result[colname], dtype='U16')
    table = Table(
        [tile_list, ra, dec, nsources,
         result['RA_min'], result['RA_max'],
         result['Dec_min'], result['Dec_max']
         ],
        names=('Tile', 'RA', 'Dec', 'NSources',
               'RA_min', 'RA_max',
               'Dec_min', 'Dec_max'
               )
    )

    """
    result['Xpos_min'], result['Xpos_max'],
    result['Ypos_min'], result['Ypos_max']
    'Xpos_min', 'Xpos_max',
    'Ypos_min', 'Ypos_max'
    """

    table.info(['stats', 'attributes'])

    dec_unique = np.unique(table['Dec'])
    print(dec_unique)


    outfile = 'tmp.fits'
    table.write(outfile, overwrite=True)

    print('Elapsed time:', time.time() - t0, 'seconds')

    xdata = table['RA']
    ydata = table['Dec']

    title = table_name
    label = str(nrows)
    plt.plot(xdata, ydata,
             marker='+', linestyle='none',
             label=label)
    plt.legend()

    plt.title(title)
    plt.xlabel('RA')
    plt.ylabel('Dec')

    saveplot=True
    plotfile = 'radec.png'
    if saveplot:
        print('Saving:', plotfile)
        plt.savefig(plotfile)

    plt.show()

    sys.exit()

    return



def tap_count_rows(table=None,
                   colnames_radec=['RAJ2000', 'DEJ2000'],
                   columnName='RAJ2000',
                   GET_COLNAMES=True,
                   GET_RADEC_LIMITS=False):
    """

    """
    import pyvo


    print()
    print(f'tap_count_rows')

    t0 = time.time()

    # Connect to VizieR TAP service
    print(f'pyvo.dal.TAPService')
    tap = pyvo.dal.TAPService("http://tapvizier.cds.unistra.fr/TAPVizieR/tap")

    # note the extra exclosing quotes needed; thank you Claude for helping me.
    # CDS currently does not support nrows info
    """
    table_info = tap.tables[f'"{table}"']
    print(f"{table}: {table_info.nrows} rows")
    """

    if GET_COLNAMES:
        query = f"SELECT TOP 1 * FROM \"{table}\""
        print(query)
        table_info = tap.search(query)
        print('table_info.to_table().colnames')
        print(table_info.to_table().colnames)


    # Execute COUNT query
    query = f'SELECT COUNT(*) AS num_rows FROM "{table}"'
    query = f'SELECT COUNT({colnames_radec[0]}) AS num_rows FROM "{table}"'
    print(query)
    result = tap.search(query)
    num_rows = result['num_rows'][0]
    print('Elapsed time:', time.time() - t0, 'seconds')
    print(f"Number of rows: {num_rows}")

    SKIP_THIS = False
    if not SKIP_THIS:
        """ This is very slow for large tables like kids_dr5 """
        t1 = time.time()
        print(f'astroquery Taplus')
        tap = TapPlus(url="http://tapvizier.cds.unistra.fr/TAPVizieR/tap")
        query = f'SELECT COUNT(*) FROM "{table}"'
        print(f'query: {query}')
        result = tap.launch_job(query).get_results()
        print(table, result)
        print('Elapsed time:', time.time() - t1, 'seconds')
        print('Total Elapsed time:', time.time() - t0, 'seconds')



        if GET_RADEC_LIMITS:
            query = f"""
            SELECT
            MIN({colnames_radec[0]}),
            MAX({colnames_radec[0]}),
            MIN({colnames_radec[1]}),
            MAX({colnames_radec[1]})
            FROM "{table}"
            """

            t1 = time.time()
            tap = pyvo.dal.TAPService("http://tapvizier.cds.unistra.fr/TAPVizieR/tap")
            print(f'query: {query}')
            result = tap.search(query)
            print(f'result: {result}')
            print('Elapsed time:', time.time() - t0, 'seconds')
            print('Elapsed time:', time.time() - t1, 'seconds')


    return


def mk_url_vizier_catalog(catalog=None):
    """
    e.g. catalog = 'IX/70'

    """

    url_base = 'https://vizier.cds.unistra.fr/viz-bin/VizieR'

    if catalog is not None:
        url = url_base + '?-source=' + catalog

    if catalog is None:
        url = url_base

    return url


def cds_vizier_xmatch(table1=None,
                      table2=None,
                      colname_xmatch_separation = 'angDist',
                      max_distance=4.0,
                      test=False):

    try_logger()

    t0 = time.time()


    cat1 = table1
    cat2 = table2

    max_distance = max_distance * u.arcsec
    print(max_distance)
    t0 = time.time()

    print(f'cat1: {cat1}')
    print(f'cat2: {cat2}')
    print(f'Maximum xmatch distance: {max_distance}')

    if test:
        table_xmatch = XMatch.query(cat1, cat2, max_distance)

    table_xmatch = XMatch.query(cat1, cat2, max_distance)
    print('Elapsed time:', time.time() - t0, 'seconds')

    print(f'Number of rows: {len(table_xmatch)}')

    table_xmatch.info()
    table_xmatch.info(['attributes', 'stats'])

    print(max_distance)
    print(f'All matchs: {len(table_xmatch)}')
    n_all = len(table_xmatch)
    print()
    itest = table_xmatch[colname_xmatch_separation] < (max_distance/2.0)
    ntest_2 = len(table_xmatch[itest])
    print(len(table_xmatch[itest]))

    itest = table_xmatch[colname_xmatch_separation] < (max_distance/4.0)
    ntest_4 = len(table_xmatch[itest])
    print(len(table_xmatch[itest]))

    itest = table_xmatch[colname_xmatch_separation] < (max_distance/8.0)
    ntest_8 = len(table_xmatch[itest])
    print(len(table_xmatch[itest]))

    print()


    # 16  - 4 = 12

    print(f'Number of matchs > {max_distance/2}: {n_all - ntest_2} {(n_all - ntest_2)/12}')
    print(f'Number of matchs > {max_distance/4}: {n_all - ntest_4}')

    print(f'Number of matchs > {max_distance/8}: {n_all - ntest_8}')

    print(f'colname_xmatch_separation: {colname_xmatch_separation}')

    return table_xmatch



def cds_vizier_xmatch_tap(table1_name=None,
                          table1_catalog=None,
                          table1_colnames_radec=['RAJ2000', 'DEJ2000'],
                          table2_name=None,
                          table2_catalog=None,
                          table2_colnames_radec=['RAJ2000', 'DEJ2000'],
                          max_distance=10.0,
                          offset_ra=None,
                          offset_dec=None,
                          get_metadata=False,
                          get_only_metadata=False,
                          async_tap=False,
                          xmatch_table_stats=True,
                          outfile_xmatch='cds_vizier_xmatch.fits',
                          outfile_xmatch_suffix=None,
                          add_nrows=True,
                          debug=False,
                          verbose=False):
    """


    """

    import time
    from datetime import datetime, timezone
    import logging


    import pyvo
    from astroquery.utils.tap.core import TapPlus

    try_logger()
    logger.info('\n')
    t0 = time.time()

    # Connect to CDS VizieR TAP service
    service = pyvo.dal.TAPService(
        "http://tapvizier.u-strasbg.fr/TAPVizieR/tap")
    if verbose:
        service.describe()
        logger.info('\n')

    #result = service.search("""SELECT table_name, description
    #    FROM TAP_SCHEMA.tables
    #    WHERE LOWER(description) LIKE '%desi%'
    #        OR LOWER(table_name) LIKE '%desi%'
    #""")
    #print(result)

    # astroquery tap
    #tap = TapPlus(url="https://TAPVizieR.cds.unistra.fr/TAPVizieR/tap")
    #help(tap)


    if DEBUG:
        logger.info('')
        input('Enter any key to continue... ')

        print([tab_name for tab_name in service.tables.keys()])
        # help(service)

    # print(service.description)

    # Get table1 metadata
    if table1_name is None:
        table1_name = 'II/383/kids_dr5'
        table2_catalog = 'II/383'

    print(f'table1_catalog: {table1_catalog}')

    if get_metadata or get_only_metadata:
        logger.info('')
        print(f'Vizier(catalog=table1_catalog).get_catalog_metadata()')
        catalog_metadata = Vizier(catalog=table1_catalog).get_catalog_metadata()
        catalog_metadata.info()
        print(catalog_metadata)

    query = f"SELECT TOP 1 * FROM \"{table1_name}\""
    print(query)
    logger.info('')
    table1_info = service.search(query)

    table1_info = table1_info.to_table()
    ncolumns1 = len(table1_info.colnames)

    print(f'table1_name: {table1_name}')
    print(f'Number of columns: {ncolumns1}')
    if verbose:
        print('table1_info.to_table().colnames')
        print(table1_info.colnames)


    if table2_name is None:
        table2_name = 'II/383/kids_dr5'
        table2_catalog = 'II/383'

    if get_metadata or get_only_metadata:
        print(f'Vizier(catalog=table2_catalog).get_catalog_metadata()')
        logger.info('')
        catalog_metadata = Vizier(catalog=table2_catalog).get_catalog_metadata()
        catalog_metadata.info()
        print(f'table2 catalog_metadata: {catalog_metadata}')

    query = f"SELECT TOP 1 * FROM \"{table2_name}\""
    print(query)
    logger.info('')
    table2_info = service.search(query)

    print(f'table2_name: {table2_name}')
    table2_info = table2_info.to_table()
    ncolumns2 = len(table2_info.colnames)
    print(f'Number of columns: {ncolumns2}')
    if verbose:
        print('table2_info.to_table().colnames')
        print(table2_info.colnames)

    if get_only_metadata:
        return

    # ADQL query to crossmatch Milliquas with KIDS DR5 within 10 arcsecond
    nrows = 100
    adql_query = f"""
    SELECT TOP {nrows}
    mq.Name AS milliquas_name,
    mq.RAJ2000 AS mq_ra,
    mq.DEJ2000 AS mq_dec,
    mq.z AS redshift,
    kids.RAJ2000 AS kids_ra,
    kids.DEJ2000 AS kids_dec,
    kids.rmagGAAP AS kids_mag,
    DISTANCE(
        POINT('ICRS', mq.{table1_colnames_radec[0]},
                      mq.{table1_colnames_radec[1]})
        POINT('ICRS', kids.RAJ2000, kids.DEJ2000)
    ) AS separation_deg
    FROM
    "{table1_name}" AS mq
    JOIN "{table2_name}" AS kids
     ON 1=CONTAINS(
          POINT('ICRS', mq.RAJ2000, mq.DEJ2000),
          CIRCLE('ICRS', kids.RAJ2000, kids.DEJ2000, 15.0/3600.0)
     )
    ORDER BY separation_deg
    """

    max_distance_deg = max_distance/3600.0

    # TOP {nrows}
    adql_query = f"""
    SELECT
    t1.*,
    t2.*,
    DISTANCE(
        POINT('ICRS', t1.{table1_colnames_radec[0]}, t1.{table1_colnames_radec[1]}),
        POINT('ICRS', t2.{table2_colnames_radec[0]}, t2.{table2_colnames_radec[1]})
    ) AS separation_deg
    FROM
    "{table1_name}" AS t1
    JOIN "{table2_name}" AS t2
     ON 1=CONTAINS(
          POINT('ICRS', t1.{table1_colnames_radec[0]}, t1.{table1_colnames_radec[1]}),
          CIRCLE('ICRS', t2.{table2_colnames_radec[0]}, t2.{table2_colnames_radec[1]}, {max_distance_deg})
     )
    ORDER BY separation_deg
    """

    if offset_ra is not None:
        """Key points:
        Dec offset: simply add offset_arcsec / 3600.0 to convert to degrees
        RA offset: divide by COS(RADIANS(dec)) to account for the spherical geometry
        — RA degrees shrink towards the poles
        If your offsets are already in degrees, drop the /3600.0.
        Note: Not all TAP services support COS and RADIANS in ADQL. If yours doesn't,
        you'd need to precompute the RA correction or apply the offset in degrees directly
        (which is only accurate for small offsets near the equator). Opus 4.5"""



        # Query the TAP capabilities

        result = service.search("SELECT * FROM TAP_SCHEMA.tables")
        print(result.to_table()['table_name'])

        # result = service.search("SELECT * FROM TAP_SCHEMA.functions")
        # print(result.to_table())

        # Offsets in arcseconds
        ra_offset_arcsec = 0
        dec_offset_arcsec = offset_dec

        ra_offset_deg = ra_offset_arcsec/3600.0
        dec_offset_deg = dec_offset_arcsec/3600.0

        adql_query = f"""
        SELECT
        t1.*,
        t2.*,
        DISTANCE(
            POINT('ICRS', t1.{table1_colnames_radec[0]},
                          t1.{table1_colnames_radec[1]} + ({dec_offset_deg})),
            POINT('ICRS', t2.{table2_colnames_radec[0]}, t2.{table2_colnames_radec[1]})
        ) AS separation_deg
        FROM
        "{table1_name}" AS t1
        JOIN "{table2_name}" AS t2
         ON 1=CONTAINS(
              POINT('ICRS', t1.{table1_colnames_radec[0]},
                            t1.{table1_colnames_radec[1]} + ({dec_offset_deg})),
              CIRCLE('ICRS', t2.{table2_colnames_radec[0]}, t2.{table2_colnames_radec[1]}, {max_distance_deg})
        )
        """


    print(f'ADQL query: {adql_query}')
    logger.info('')
    # Execute the query
    print("Running TAP crossmatch query...")

    if hasattr(service, 'timeout'):
        print(f"Service timeout: {service.timeout}")
    if hasattr(service, 'hardlimit'):
        print(f"Service hardlimit: {service.hardlimit}")
    if hasattr(service, 'maxrec'):
        print(f"Service maxrec: {service.maxrec}")

    print(f'table1_name: {table1_name}')
    print(f'table2_name: {table2_name}')
    print(f'max_distance (arcsec): {max_distance}')
    print('Elapsed time:', time.time() - t0, 'seconds')
    logger.info('')

    print(service.describe)

    logger.info('')
    print('Elapsed time:', time.time() - t0, 'seconds')
    logger.info(f'async TAP query: {async_tap}')

    if async_tap:

        t0_tap = time.time()
        #job = service.run_async(adql_}query)
        job = service.submit_job(adql_query)
        logger.info(f'Async TAP job submitted: {job.url}')
        logger.info(f'Async TAP job phase: {job.phase}')
        job.run()

        # Wait for completion (polls automatically)
        job.wait()
        logger.info(f'Async TAP Job finished; phase: {job.phase}')
        logger.info(f'Job elapsed time: {time.time() - t0_tap} seconds')
        print('Elapsed time:', time.time() - t0, 'seconds')

        explore_result = False
        if explore_result:
        # help(job)
            # Check result size
            # print(f"Result size: {job.result_size}")  # bytes, if available
            # Get just the column info without downloading all data
            logger.info(f"Result URI: {job.result_uri}")

            import requests

            response = requests.get(job.url)
            print(f'job.url: {response.text}')

            # Stream just the first part to get column definitions
            response = requests.get(job.result_uri, stream=True)
            header_chunk = b''
            for chunk in response.iter_content(chunk_size=4096):
                header_chunk += chunk
                if b'<TABLEDATA>' in header_chunk or b'<BINARY>' in header_chunk:
                    break
                if len(header_chunk) > 50000:  # Safety limit
                    break
            response.close()

            # Parse just the header
            header_text = header_chunk.decode('utf-8', errors='ignore')

            # Find column definitions
            import re
            fields = re.findall(r'<FIELD[^>]*name="([^"]*)"[^>]*datatype="([^"]*)"', header_text)
            print(f"Columns ({len(fields)}):")
            for name, dtype in fields:
                print(f"  {name}: {dtype}")

            logger.info(f'Job elapsed time: {time.time() - t0_tap} seconds')
            response = requests.head(job.result_uri)
            print(f"Content-Length: {response.headers.get('Content-Length')} bytes")
            print(f"Content-Type: {response.headers.get('Content-Type')}")
            logger.info(f'Job elapsed time: {time.time() - t0_tap} seconds')

            import requests

            headers = {'Range': 'bytes=0-50000'}
            response = requests.get(job.result_uri, headers=headers)
            header_text = response.text

            logger.info(f'Job elapsed time: {time.time() - t0_tap} seconds')
            input('Enter any key to continue... ')

            print('header_text[:5000]')
            #print(header_text[:20000])  # Inspect the structure
            help(header_text)
            print(len(header_text))
            logger.info('')
            input('Enter any key to continue... ')

            print(header_text)  # Inspect the structure

        # Get results
        logger.info(f'Fetch results')
        results = job.fetch_result()
        print(f"Results size: {sys.getsizeof(results)} bytes")


    if not async_tap:
        #help(service.search)
        results = service.search(adql_query)

    logger.info('')
    print('Elapsed time:', time.time() - t0, 'seconds')
    nrows = len(results)
    print(f'Query finished {nrows} rows')
    if debug:
        logger.info('')
        input('Enter any key to continue... ')


    # Convert to astropy table
    table = results.to_table()
    print(f'Query results converted to table')
    print('Elapsed time:', time.time() - t0, 'seconds')

    # Optionally save to file as vot
    #table.write('milliquas_kids_crossmatch.vot', format='votable',
    #            overwrite=True)

    if verbose:
        table.info()

    RUN_table_info_stats = False
    if RUN_table_info_stats:
        table.info(['attributes', 'stats'])
        table.info()

    if verbose:
        print("\ncolumn dtypes:")
        print(table.dtype)

    logger.info('')

    # Find all object/string columns
    object_columns = [col for col in table.colnames if table[col].dtype == 'object']

    print("Object-type columns and their maximum lengths:")
    print("-" * 50)

    # Compute and display maximum length for each object column
    for col in object_columns:
        # Filter out None/NaN values and compute max length
        non_null_values = [str(val) for val in table[col] if val is not None]
        if non_null_values:
            max_len = max(len(val) for val in non_null_values)
            print(f"{col}: max length = {max_len}")

            # Convert to fixed-length string dtype
            table[col] = table[col].astype(f'U{max_len}')
        else:
            print(f"{col}: all values are None/null")

    if verbose:
        print("\nUpdated dtypes:")
        print(table.dtype)

    if verbose:
        table.info()
        logger.info('')
        input('Enter any key to continue... ')

    if xmatch_table_stats:
        table.info(['attributes', 'stats'])
        table[table1_colnames_radec[0]].info(['attributes', 'stats'])
        table[table1_colnames_radec[1]].info(['attributes', 'stats'])
        table.info()

    # Display results
    print(f"\nFound {len(table)} matches within {max_distance} arcseconds")
    if verbose:
        table.pprint(max_lines=20, max_width=120)
        table.info()

    table_fix_column_units(table=table)

    # Optionally save to file
    #filename_add_nrows = True
    #if filename_add_nrows:

    path = Path(outfile_xmatch)
    if outfile_xmatch_suffix is  None:
        outfile = f"{path.stem}_{nrows}{path.suffix}"
    if outfile_xmatch_suffix is not None:
        outfile = f"{path.stem}_{outfile_xmatch_suffix}_{nrows}{path.suffix}"

    print(f'outfile: {outfile}')

    now = datetime.now(timezone.utc)
    # 2026-02-28 13:45:00.1 UTC
    print(now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-5] + " UTC")
    print(f'Elapsed time: {(time.time() - t0):.3f} seconds; {now.strftime("%Y-%m-%d %H:%M:%S.%f")} UTC')

    t1 = time.time()
    logger.info(f'Writing {outfile}')
    table.write(outfile, overwrite=True)
    logger.info(f'{outfile} written successfully')
    print('Table write elapsed time:', time.time() - t1, 'seconds')
    print('Total elapsed time:', time.time() - t0, 'seconds')

    print()
    print(f'table1_name: {table1_name}')
    print(f'table2_name: {table2_name}')
    print(f'{len(table)} rows')
    print(f'max_distance (arcsec): {max_distance}')
    logger.info('')

    table_test = table
    print()
    print(f'{table1_name}')
    print(table1_colnames_radec)
    print(table[table1_colnames_radec[0]].unit)
    print(table[table1_colnames_radec[1]].unit)

    print(f'RA min: {np.min(table_test[table1_colnames_radec[0]])}')
    print(f'RA max: {np.max(table_test[table1_colnames_radec[0]])}')
    print(f'Dec min: {np.min(table_test[table1_colnames_radec[1]])}')
    print(f'Dec max: {np.max(table_test[table1_colnames_radec[1]])}')

    print()
    print(f'{table2_name}')
    print(table2_colnames_radec)
    print(table[table2_colnames_radec[0]].unit)
    print(table[table2_colnames_radec[1]].unit)
    print(f'RA min: {np.min(table_test[table2_colnames_radec[0]])}')
    print(f'RA max: {np.max(table_test[table2_colnames_radec[0]])}')
    print(f'Dec min: {np.min(table_test[table2_colnames_radec[1]])}')
    print(f'Dec max: {np.max(table_test[table2_colnames_radec[1]])}')

    logger.info('')

    test_source_radec1 = table[table1_colnames_radec][0]
    test_source_radec2 = table[table2_colnames_radec][0]
    print(f'test_source_radec1: {test_source_radec1}')
    print(f'test_source_radec2: {test_source_radec2}')

    table1_test = tap_table_check(table=None, radec=None)
    table2_test = tap_table_check(table=None, radec=None)

    logger.info('')
    input('Enter any key to continue... ')

    print(f'table: {len(table)}')
    colname_xmatch_separation = 'separation_deg'
    itest = table[colname_xmatch_separation] < (max_distance/(2.0*3600.0))
    ntest = len(table[itest])
    print(f'table: {len(table[itest])}')

    itest = table[colname_xmatch_separation] < (max_distance/(4.0*3600.0))
    ntest = len(table[itest])
    print(len(table[itest]))

    use_fitsio = False
    if use_fitsio:
        # try fitsio option
        t1 = time.time()
        fitsio.write(outputfile_xmatch, table.as_array(), clobber=True)
        table.write(fitsio.write(outputfile_xmatch, clobber=True))
        print('Table write elapsed time:', time.time() - t1, 'seconds')
        print('Total elapsed time:', time.time() - t0, 'seconds')
    print()

    # Read back test
    read_back_test = False
    if read_back_test:
        table_test = Table.read(outfile)
        table_test.info(['attributes', 'stats'])

    result_radec_limits(table=table_test,
                        table1_name=table1_name,
                        table2_name=table2_name,
                        table1_colnames_radec=table1_colnames_radec,
                        table2_colnames_radec=table2_colnames_radec
                        )

    logger.info('')
    t0 = time.time()

    print(f'astroquery Taplus')
    input('Enter any key to continue... ')

    # There is a limit of 2000 rows. If you need more than that, you must use asynchronous queries.
    tap = TapPlus(url="http://tapvizier.cds.unistra.fr/TAPVizieR/tap")
    print(f'query: {adql_query}')
    result = tap.launch_job(adql_query).get_results()
    help(result)
    print(f'{len(result)} rows')
    print(result)

    logger.info('')
    print('Elapsed time:', time.time() - t0, 'seconds')

    table_test = result
    print()
    print(f'{table1_name}')
    print(table1_colnames_radec)
    print(table[table1_colnames_radec[0]].unit)
    print(table[table1_colnames_radec[1]].unit)

    print(f'RA min: {np.min(table_test[table1_colnames_radec[0]])}')
    print(f'RA max: {np.max(table_test[table1_colnames_radec[0]])}')
    print(f'Dec min: {np.min(table_test[table1_colnames_radec[1]])}')
    print(f'Dec max: {np.max(table_test[table1_colnames_radec[1]])}')

    print()
    print(f'{table2_name}')
    print(table2_colnames_radec)
    print(table[table2_colnames_radec[0]].unit)
    print(table[table2_colnames_radec[1]].unit)
    print(f'RA min: {np.min(table_test[table2_colnames_radec[0]])}')
    print(f'RA max: {np.max(table_test[table2_colnames_radec[0]])}')
    print(f'Dec min: {np.min(table_test[table2_colnames_radec[1]])}')
    print(f'Dec max: {np.max(table_test[table2_colnames_radec[1]])}')

    logger.info('')

    return table


def tap_table_check(catalog="VII/289",
                    table="VII/289/dr16q",
                    table_colnames_radec=['RAJ2000', 'DEJ2000'],
                    radec_centre=[180.0, -60.0],
                    box_width_deg=(1.0/60),
                    test=False,
                    verbose=False):
    """
    exmaple:

    table_test = tap_table_check(table=table1_name,
                                 table_colnames_radec=table1_colnames_radec,
                                 radec_centre=[180.0, -60.0])


    See `the vizier constraints page
 |      https://vizier.cds.unistra.fr/vizier/vizHelp/syntax.htx  for details.

    """
    import sys
    import time
    from astroquery.vizier import Vizier
    help(Vizier)

    t0 = time.time()

    logger = try_logger()
    logger.info('')

    metadata = Vizier(catalog=catalog).get_catalog_metadata()
    print(type(metadata))
    print(len(metadata))
    print()
    print(metadata)

    print()


    # Vizier interface

    test = True
    if test:

        print('See https://vizier.cds.unistra.fr/vizier/vizHelp/syntax.htx')
        # need to be careful about spherical distortions
        radec_centre = [180.0, 0.0]
        box_width_deg = (0.5)

        dec_min = radec_centre[1] - (box_width_deg/2.0)
        dec_max = radec_centre[1] + (box_width_deg/2.0)

        dec_abs_max = max(abs(dec_min), abs(dec_max))
        print(dec_abs_max, np.deg2rad(dec_abs_max))
        cosdec = np.cos(np.deg2rad(dec_abs_max))

        print(f'dec_abs_dec, cosdec, 1/cosdec: '
              f'{dec_abs_max}, {cosdec}, {1.0/cosdec}')

        ra_min_deg = radec_centre[0] - (box_width_deg/(2.0*cosdec))
        ra_max_deg = radec_centre[0] + (box_width_deg/(2.0*cosdec))

        dec_min_deg = radec_centre[1] - (box_width_deg/(2.0))
        dec_max_deg = radec_centre[1] + (box_width_deg/(2.0))

        print(f'RA limits: {ra_min_deg}, {ra_max_deg}')
        print(f'Dec limits: {dec_min_deg}, {dec_max_deg}')

        # sys.exit()

        # RAJ2000=">=179.9 & <=180.1"
        RAJ2000_Constraint = f'>= {ra_min_deg} & <= {ra_max_deg}'
        # RAJ2000_Constraint = ">= 179.9 & <= 180.1"

        #DEJ2000=">=-0.1 & <= 0.1",
        DEJ2000_Constraint = f'>= {dec_min_deg} & <= {dec_max_deg}'
        # DEJ2000_Constraint = ">= -0.1 & <= 0.1"

        print(f'RAJ2000_Constraint: {RAJ2000_Constraint} '
              f'{type(RAJ2000_Constraint)} '
              f'{len(RAJ2000_Constraint)}')
        print(f'DEJ2000_Constraint: {DEJ2000_Constraint} '
              f'{type(DEJ2000_Constraint)} '
              f'{len(DEJ2000_Constraint)}')

        vizier = Vizier(row_limit=99)
        vizier = Vizier(columns=["**"])
        # Just fetch one row to inspect the table structure
        catalog = "VII/289"
        table = "VII/289/dr16q"
        result_vizier = vizier.query_constraints(
            catalog=table,
            RAJ2000=RAJ2000_Constraint,
            #RAJ2000=">=179.9 & <= 180.1",
            DEJ2000=DEJ2000_Constraint,
            #DEJ2000=">=-0.1 & <= 0.1",
            row_limit=99)

        # sys.exit()

    # return a Table list for the catalog which is confusing;
    # row limit applies to the list of tables!
    logger.info('')
    print(f'Result type {type(result_vizier)}')
    print('Elapsed time:', time.time() - t0, 'seconds')

    print(f'len(result_vizier): {len(result_vizier)}')

    if len(result_vizier) > 1:
        logger.info('To many tables returned; exiting...')
        sys.exit()

    if len(result_vizier) == 1:
        result_vizier = result_vizier[0]
    logger.info('')
    print(result_vizier.colnames)
    print(f'Number of columns: {len(result_vizier.colnames)}')
    logger.info('')

    result_vizier.info('attributes')
    logger.info('')
    result_vizier.info('stats')

    logger.info('')
    print(f'Number of rows: {len(result_vizier)}')
    print(f'Number of columns: {len(result_vizier.colnames)}')
    print(result_vizier.colnames)
    logger.info('')
    input('Enter any key to continue... ')

    t0 = time.time()
    # Connect to VizieR TAP service
    print(f'Use pyvo.dal.TAPService')
    tap_service = \
        pyvo.dal.TAPService("http://tapvizier.cds.unistra.fr/TAPVizieR/tap")

    #tap_service.describe()

    # note the extra exclosing quotes needed; thank you Claude for helping me.
    # CDS currently does not support nrows info
    """
    table_info = tap.tables[f'"{table}"']
    print(f"{table}: {table_info.nrows} rows")
    """

    query = f"SELECT TOP 1 * FROM \"{table}\""
    print(f'query:\ {query}')
    tap_result = tap_service.search(query)
    print('table_result.to_table().colnames')
    print(tap_result.to_table().colnames)
    print(f'Number of columns: {len(tap_result.to_table().colnames)}')
    logger.info('')

    query = f"""
        SELECT *
        FROM \"{table}\"
        WHERE
            {table_colnames_radec[0]} >= {ra_min_deg} AND
            {table_colnames_radec[0]} <= {ra_max_deg} AND
            {table_colnames_radec[1]} >= {dec_min_deg} AND
            {table_colnames_radec[1]} <= {dec_max_deg}
        """
    print(f'query:\ {query}')
    result = tap_service.search(query)

    print(f'type(result): {type(result)}')
    print(f'Number of rows: {len(result)}')
    print('result.to_table().colnames')
    print(result.to_table().colnames)
    print(f'Number of columns: {len(result.to_table().colnames)}')
    logger.info('')
    print('Elapsed time:', time.time() - t0, 'seconds')

    return result


def fix_votable_object(table, dtype_new='bool',
                       verbose=False):
    """
    convert the columns with dtype = object which is not
    supported by FITs to bool
    """

    from astropy.table import Table, Column

    for (icol, column) in enumerate(table.columns):
        if verbose:
            print(icol, table.columns[icol].name,
                table.columns[icol].format,
                table.columns[icol].dtype)

        if table.columns[icol].dtype == 'object':
            colname = table.columns[icol].name
            NewColumn = Table.Column(table[colname].data, dtype='bool')
            table.replace_column(colname, NewColumn)

    return table


def table_fix_column_units(table=None):
    """
    fix the bad column units in the table
    """
    for icol, col in enumerate(table.colnames):
        try:
            # Test if this column's unit is FITS-compatible
            from astropy.io.fits import column
            if hasattr(table[col], 'unit') and table[col].unit is not None:
                str(table[col].unit.to_string('fits'))
        except Exception:
            print(f"Removing non-standard unit from column: {icol} {col}")
            print(table[col].unit)
            table[col].unit = None

    return table



def count_rows_table(table=None, GET_COLNAMES=True):
    """

    """
    if GET_COLNAMES:
        query = f"SELECT TOP 1 * FROM \"{table}\""
        print(query)
        table_info = service.search(query)
        print('table_info.to_table().colnames')
        print(table_info.to_table().colnames)


    t1 = time.time()
    print('Elapsed time:', time.time() - t0, 'seconds')
    tap_count_rows(table=table_name,
                   colnames_radec=table_colnames_radec)
    print('Elapsed time:', time.time() - t1, 'seconds')


    return


def get_table_lists(DEBUG=False):
    """
    need to define concepts of Vizier catalog and table

    Could also get the bibcode from CDS or print the CDS page
    See: https://cdsarc.cds.unistra.fr/viz-bin/cat/IX/70
    #
    #
    # 0: VII/294/catalog',
    # The Million Quasars (Milliquas) catalogue, version 8 (Flesch, 2023)
    #
    # 1: 'IX/70/csc21mas'; Chandra Source Catalog Release 2 (CSC 2.1) (Evans+, 2024)
    # https://tapvizier.cds.unistra.fr/adql/?IX/70
    # http://vizier.cds.unistra.fr/viz-bin/VizieR-?-source=IX/70/
    # TAP: RA_ICRS, DE_ICRS

    # 2: IX/69/xmm4d13s','The 4XMM-DR13 Catalog, slim version'; Standard/Original',
    # https://tapvizier.cds.unistra.fr/adql/?IX/69
    # TAP: RA_ICRS, DE_ICRS

    #'J/ApJS/255/30/comp',
    #'J/A+A/661/A1/main',
    #'J/A+A/661/A3/ctpmain']

    # The Sloan Digital Sky Survey Quasar Catalog Sixteenth Data Release (2020)
    # https://cdsarc.cds.unistra.fr/viz-bin/cat/VII/289
    """


    # Vizier table name
    table_name_list = ['VII/294/catalog',   # 0 MilliquasV8
                       'IX/70/csc21mas',    # 1 CSC21
                       'IX/69/xmm4d13s',    # 2 4XMMDR13
                       'IX/66/xmm411st',    # 3 4XMMDR11st
                       'J/ApJS/255/30/comp',# 4 VLASS_QL1
                       'J/A+A/661/A1/main', # 5
                       'J/A+A/661/A3/ctpmain',# 6
                       'VII/289/dr16q',     # 7
                       'I/356/qsocand',     # 8
                       'II/383/kids_dr5',   # 9 kids_dr5
                       'IX/71/xmmsl3c',
                       'V/161/zcatdr1', #11
                       'I/359/vhs_dr4', # 12 VHS dR4
                       'I/355/gaiadr3'] # 13 Gaia DR3

    # by hand from count queries
    table_count_rows = [1021800, # 0
                        407806,  # 1
                        656997,  # 2
                        1780775, # 3
                        3381277, # 4
                        27910,   # 5
                        27369,   # 6
                        750414,  # 7
                        6649162, # 8
                        138812117, # 9
                        -99,# 10
                        -99,# 11
                        -99, # 12
                        -99] # 13

    table_count_columns = [17, # 0
                        486, # 1
                        -99, # 2
                        -99, # 3
                        -99, # 4
                        -99, # 5
                        -99, # 6
                        -99, # 7
                        -99, # 8
                        -99, # 9
                        -99, # 10
                        -99, # 11
                        -99, # 12
                        -99] # 13

    # Vizier catalog name; could be derived from the table name
    table_catalog_list =['VII/294', # 0
                         'IX/70',
                         'IX/69',
                         'IX/66',
                         'J/ApJS/255/30',
                         'J/A+A/661/A1',
                         'J/A+A/661/A3',
                         'VII/289',
                         'I/356', #8
                         'II/383',   # 9
                         'IX/71', # 10
                         'V/161', #11
                         'I/359', # 12
                         'I/355'] # 13



    table_label_list = ['MQv8_Flesch+2023',
                        'csc21_Evans+2024',
                        '4xmmdr13s_Webb+2023',
                        '4xmmdr11st_Traulsen+2022',
                        'VLASS_QL1_Gordon+2021',
                        'eFEDS_Brunner+2022',
                        'eFEDS_Salvato+2022',
                        'sdss_dr16q_Lyke+2020',
                        'gaiadr3_qsocand_Gaia+2022',
                        'kids_dr5_Wright+2024',
                        'XMMSL3_XMM-SSC+2025',
                        'DESI_DR1_zcat+2025', #11
                        'VHS_DR4+2019', # 12 VHS dR4
                        'gaiadr3'] # 13



    table_colnames_radec = [['RAJ2000', 'DEJ2000'], # 0
                            ['RAICRS', 'DEICRS'],
                            ['RA_ICRS', 'DE_ICRS'],
                            ['RAJ2000', 'DEJ2000'],
                            ['RAJ2000', 'DEJ2000'], # 4
                            ['RA_ICRS', 'DE_ICRS'],
                            ['RAc', 'DEc'],
                            ['RAJ2000', 'DEJ2000'],
                            ['RA_ICRS', 'DE_ICRS'], # 8
                            ['RAJ2000', 'DEJ2000'],
                            ['RAJ2000', 'DEJ2000'], # 10
                            ['RAICRS', 'DEICRS'], # 11
                            ['RAJ2000', 'DEJ2000'], #12
                            ['RAJ2000', 'DEJ2000']] # 13


    table_metadata = [table_count_rows, table_count_columns]

    return table_name_list, table_catalog_list, table_label_list, \
        table_colnames_radec, table_metadata



def set_xmatch_input_tables(itable1=0,
                            itable2=1,
                            LIST_TABLES_ONLY=False):

    try_logger()
    logger.info('')

    table_name_list, table_catalog_list, table_label_list, table_colnames_radec_list, \
        table_metadata = get_table_lists()

    table_count_rows = table_metadata[0]
    table_count_columns = table_metadata[1]

    table_name = ['name1', 'name2']
    table_catalog = ['catalog1', 'catalog2']
    table_label = [table_label_list[0], table_label_list[1]]

    """
    table_name[1] = 'II/383/kids_dr5'
    table_catalog[1] = 'II/383'

    table_name[1] = 'VII/289/dr16q'
    table_catalog[1] = 'VII/289'
    """


    print(f"Available tables:")
    print(f'table_name_list: {table_name_list}')
    for itable, table in enumerate(table_name_list):
        print(f"{itable}: "
              f"{table_label_list[itable]}; "
              f"{table_name_list[itable]}; "
              f"{table_catalog_list[itable]}; "
              f"{table_colnames_radec_list[itable]}; "
              f"{table_count_rows[itable]}; "
              f"{table_count_columns[itable]}")

    logger.info('')

    if LIST_TABLES_ONLY:
        return

    table_catalog[0] = table_catalog_list[itable1]
    table_label[0] = table_label_list[itable1]
    table_name[0] = table_name_list[itable1]


    table_catalog[1] = table_catalog_list[itable2]
    table_label[1] = table_label_list[itable2]
    table_name[1] = table_name_list[itable2]


    table_colnames_radec = [table_colnames_radec_list[itable1], table_colnames_radec_list[itable2]]

    print(f'itable1: {itable1}')
    print(f'itable2: {itable2}')
    print(f'table_name: {table_name}')
    print(f'table_catalog: {table_catalog}')
    print(f'table_label: {table_label}')
    print(f'table_colnames_radec: {table_colnames_radec}')


    return table_name, table_catalog, table_label, table_colnames_radec, table_metadata
        # table_count_rows, table_count_columns


def result_radec_limits(table=None,
                        table1_name=None,
                        table2_name=None,
                        table1_colnames_radec=['RAJ2000', 'DEJ2000'],
                        table2_colnames_radec=['RAJ2000', 'DEJ2000']):

    try_logger()
    logger.info('')

    # Print all metadata for a specific column
    column_name_list = table1_colnames_radec + table1_colnames_radec
    for column_name in column_name_list:
        print(table[column_name].info)
        table[column_name].info

    logger.info('')

    # Print all table metadata
    print(f'table.meta: {table.meta}')
    # It's a dictionary, so you can loop through:
    for key, value in table.meta.items():
        print(f"{key}: {value}")

    print()

    print(f'{table1_name}')
    print(table1_colnames_radec)
    print(table[table1_colnames_radec[0]].unit)
    print(table[table1_colnames_radec[1]].unit)
    print(f'RA min: {np.min(table[table1_colnames_radec[0]])}')
    print(f'RA max: {np.max(table[table1_colnames_radec[0]])}')
    print(f'Dec min: {np.min(table[table1_colnames_radec[1]])}')
    print(f'Dec max: {np.max(table[table1_colnames_radec[1]])}')

    print()
    print(f'{table2_name}')
    print(table2_colnames_radec)
    print(table[table2_colnames_radec[0]].unit)
    print(table[table2_colnames_radec[1]].unit)
    print(f'RA min: {np.min(table[table2_colnames_radec[0]])}')
    print(f'RA max: {np.max(table[table2_colnames_radec[0]])}')
    print(f'Dec min: {np.min(table[table2_colnames_radec[1]])}')
    print(f'Dec max: {np.max(table[table2_colnames_radec[1]])}')

    return


def parse_arguments():

    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run CDS Vizier TAP xmatch for two tables; \n note in CDS Vizier a catalog consists of multiple tables with table names "
    )

    parser.add_argument(
        "--list_tables",
        action="store_true",
        default=False,
        help="List available tables and exit (default: %(default)s)",
    )

    parser.add_argument(
        "--count_rows_tables",
        action="store_true",
        default=False,
        help="Count rows in all tables and exit (default: %(default)s)",
    )


    parser.add_argument(
        "--get_metadata",
        action="store_true",
        default=False,
        help="Get the table metadata (default: %(default)s)",
    )

    parser.add_argument(
        "--max_distance",
        type=float,
        default=2.0,
        help=f"xmatch maximum radius in arc seconds (default: %(default)s)",
    )

    parser.add_argument(
        "--offset_distance",
        type=float,
        default=None,
        help=f"offset distance in arc seconds (default: %(default)s)",
    )

    parser.add_argument(
        "--itable1",
        type=int,
        default=0,
        help=f"Index of the table 1 in query (default: %(default)s)",
    )

    parser.add_argument(
        "--itable2",
        type=int,
        default=1,
        help=f"Index  of table 2 in query (default: %(default)s)",
    )

    parser.add_argument(
        "--async_tap",
        action="store_true",
        default=False,
        help="Use Vizier TAP async service (default: %(default)s)",
    )

    parser.add_argument(
        "--xmatch_cds",
        action="store_true",
        default=False,
        help="Use Vizier xmatch service (default: %(default)s)",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Debug option (default: %(default)s)",
    )


    parser.add_argument(
        "--stage_debug",
        action="store_true",
        default=False,
        help="Stage debug option (default: %(default)s)",
    )


    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose option (default: %(default)s)",
    )



    parser.add_argument(
        "--pause",
        action="store_true",
        default=False,
        help="Pause at various points (default: %(default)s)",
    )


    return parser.parse_args()




def main(DEBUG=False, tap_catalog=None, RUN_COUNT_ROWS=False):
    """


    """


    from datetime import datetime, timezone
    import sys

    logger = get_logger()

    logger.info('')
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Script: {sys.argv[0]}")
    logger.info(f"Arguments: {sys.argv[1:]}")

    t0 = time.time()

    now = datetime.now(timezone.utc)
    # 2026-02-28 13:45:00.1 UTC
    print(now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-5] + " UTC")
    print(f'Elapsed time: {(time.time() - t0):.3f} seconds; {now.strftime("%Y-%m-%d %H:%M:%S.%f")} UTC')

    print(f'Elapsed time: {(time.time() - t0):.3f} seconds; {now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]} UTC')

    print(
    f'Elapsed time: {(time.time() - t0):.3f} seconds; '
    f'{now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]} UTC')

    print()
    args = parse_arguments()
    for name, value in vars(args).items():
        print(f'CLI arg: {name} = {value}')
    print()

    debug = args.debug
    xmatch_cds = args.xmatch_cds
    async_tap = args.async_tap

    max_distance = args.max_distance
    print(f'max_distance: {max_distance}')

    offset_ra = args.offset_distance
    print(f'offset_ra: {offset_ra}')

    offset_dec = args.offset_distance
    print(f'offset_dec: {offset_dec}')

    itable1 = getattr(args, 'itable1')
    print(f'itable1: {itable1}')

    itable2 = getattr(args, 'itable2')
    print(f'itable2: {itable2}')

    debug = args.debug
    pause = args.pause
    verbose = args.verbose

    GET_METADATA = args.get_metadata
    LIST_TABLES = getattr(args, 'list_tables')

    print(f'LIST_TABLES: {LIST_TABLES}')
    if LIST_TABLES:
        set_xmatch_input_tables(LIST_TABLES_ONLY=LIST_TABLES)
        sys.exit()


    COUNT_ROWS_TABLES= getattr(args, 'count_rows_tables')
    if COUNT_ROWS_TABLES:
        t1 = time.time()
        table_name_list, table_label_list, table_catalog_list, table_colnames_radec_list = \
            get_table_lists()

        for itable, table_name in enumerate(table_name_list):
            table_colnames_radec = table_colnames_radec_list[itable]
            print()
            print(itable, table_name, table_colnames_radec)
            print('Elapsed time:', time.time() - t0, 'seconds')
            tap_count_rows(table=table_name,
                       colnames_radec=table_colnames_radec)
            print('Elapsed time:', time.time() - t1, 'seconds')

        sys.exit()


    XMATCH_GET_AVAILABLE_TABLES = False
    if XMATCH_GET_AVAILABLE_TABLES:
        table_list = XMatch.get_available_tables()
        print(f'Number of CDS xmatch tables: {len(table_list)}')
        logger.info('')
        input('Enter any key to continue... ')

        tablename_substring = 'gaia'
        tablename_substring = 'desi'
        tablename_substring = 'V/161/'
        #tablename_substring = 'zcat'
        print(f'Find names with table name substring: {tablename_substring}')
        indices = [i for i, item in enumerate(table_list) if
                   tablename_substring in item.lower()]
        print(len(indices))
        print(indices)

        logger.info('')
        input('Enter any key to continue... ')

        table_list = [table_list[i] for i in indices]
        print(table_list)

        logger.info('')
        input('Enter any key to continue... ')

    #help(table_list)
    #print(table_list)


    if DEBUG:
        help(Vizier)
        help(Vizier.get_catalog_metadata)
        # cat1 = cat_gaia_dr3
        #cat1 = 'VII/294'
        #print(f' Get catalog metadata: {cat1}')
        #catalog_metadata = Vizier(catalog=cat1).get_catalog_metadata()
        #if DEBUG:
        #    help(catalog_metadata)
        # catalog_metadata.info()

    outfile_xmatch_suffix = f"an{int(max_distance)}as"

    print(f'itable1: {itable1}')
    print(f'itable2: {itable2}')
    table_names, table_catalogs, table_labels, tables_colnames_radec, \
        table_metadata = \
        set_xmatch_input_tables(itable1=itable1, itable2=itable2)

    catalog1_url = mk_url_vizier_catalog(catalog=table_catalogs[0])
    logger.info(f'catalog1 url: {catalog1_url}')
    catalog2_url = mk_url_vizier_catalog(catalog=table_catalogs[1])
    logger.info(f'catalog2 url: {catalog2_url}')

    #table_count


    table1_label = table_labels[0]
    table2_label = table_labels[1]

    outfile_xmatch = f"{table1_label}_xmatch_tap_{table2_label}.fits"

    table1_name = table_names[0]
    table2_name = table_names[1]

    table1_catalog = table_catalogs[0]
    table2_catalog = table_catalogs[1]

    table1_colnames_radec = tables_colnames_radec[0]
    table2_colnames_radec = tables_colnames_radec[1]

    logger.info('')
    DO_TABLE_CHECK = False
    if DO_TABLE_CHECK:
        table_test = tap_table_check(
            table=table1_name,
            table_colnames_radec=table1_colnames_radec,
            radec_centre=[180.0, 0.0])

        input('Enter any key to continue... ')
        sys.exit()


    if debug or verbose:
        catalog1_ivoid = f"ivo://CDS.VizieR/{table1_catalog}"
        voresource1 = pyvo.registry.search(ivoid=catalog1_ivoid)[0]
        voresource1.describe(verbose=True)
        tables = voresource1.get_tables()
        print(f"In this catalogue, we have {len(tables)} tables.")
        for table_name, table in tables.items():
            print(f"{table_name}: {table.description}")

        logger.info('')
        catalog2_ivoid = f"ivo://CDS.VizieR/{table2_catalog}"
        voresource2 = pyvo.registry.search(ivoid=catalog2_ivoid)[0]
        voresource2.describe(verbose=True)
        tables = voresource2.get_tables()
        print(f"In this catalogue, we have {len(tables)} tables.")
        for table_name, table in tables.items():
            print(f"{table_name}: {table.description}")
        logger.info('')

    if verbose:
        logger.info('')
        input('Enter any key to continue... ')


    if xmatch_cds:
        colname_xmatch_separation = 'angDist'
        table_xmatch = cds_vizier_xmatch(
            table1=table1_name,
            table2=table2_name,
            max_distance=max_distance,
            colname_xmatch_separation=colname_xmatch_separation)
        print('Elapsed time:', time.time() - t0, 'seconds')
        table_xmatch.info()
        print()
        nrows_xmatch = len(table_xmatch)
        print(f'Number of rows in xmatch result table: {nrows_xmatch}')
        ncolumns_xmatch = len(table_xmatch.colnames)
        print(f'Number of columns in xmatch result table: {ncolumns_xmatch}')
        print()
        print('Elapsed time:', time.time() - t0, 'seconds')

        table_fix_column_units(table=table_xmatch)

        outfile = 'tmp_xmatch.fits'
        print(f'Write: {outfile}')
        table_xmatch.write(outfile, overwrite=True)
        print(f'Number of rows: {len(table_xmatch)}')
        print('Elapsed time:', time.time() - t0, 'seconds')

        logger.info('')


        sys.exit()



    if GET_METADATA:
        print(f'Running cds_vizier_xmatch_tap get_only_metadata')
        cds_vizier_xmatch_tap(table1_name=table1_name,
                              table1_catalog=table1_catalog,
                              table2_name=table2_name,
                              table2_catalog=table2_catalog,
                              get_only_metadata=True,
                              debug=debug)

        sys.exit()


    print(f'RUN_COUNT_ROWS: {RUN_COUNT_ROWS}')
    if RUN_COUNT_ROWS:
        t1 = time.time()
        print('Elapsed time:', time.time() - t0, 'seconds')
        tap_count_rows(table=table1_name,
                       colnames_radec=table1_colnames_radec)
        print('Elapsed time:', time.time() - t1, 'seconds')

        t1 = time.time()
        print('Elapsed time:', time.time() - t0, 'seconds')
        tap_count_rows(table=table2_name,
                       colnames_radec=table2_colnames_radec)
        print('Elapsed time:', time.time() - t1, 'seconds')

        # sys.exit()


    cds_vizier_xmatch_tap(table1_name=table1_name,
                          table1_catalog=table1_catalog,
                          table1_colnames_radec=table1_colnames_radec,
                          table2_name=table2_name,
                          table2_catalog=table2_catalog,
                          table2_colnames_radec=table2_colnames_radec,
                          max_distance=max_distance,
                          async_tap=async_tap,
                          offset_ra=offset_ra,
                          offset_dec=0.0,
                          outfile_xmatch=outfile_xmatch,
                          outfile_xmatch_suffix=outfile_xmatch_suffix,
                          get_only_metadata=False,
                          debug=debug)

    logger.info('')
    sys.exit()

    # Case-insensitive search
    # Get indices of matching elements
    tablename_substring = 'gaia'
    print(f'Find names with table name substring: {tablename_substring}')
    indices = [i for i, item in enumerate(table_list) if
    tablename_substring in item.lower()]
    print(len(indices))
    print(indices)

    table_list = [table_list[i] for i in indices]
    print(table_list)

    cat_kids_dr5 = 'II/383/kids_dr5'
    cat_gaia_dr3 = 'I/355/gaiadr3'

    # Milliquas
    # The Million Quasars (Milliquas) catalogue, version 8 (Flesch, 2023)
    cat_MQ = 'VII/294/catalog'

    # XMM
    cat_4xmmdr10s = 'IX/63/xmm4d10s'
    cat_4xmmdr10st = 'IX/64/xmm410st'
    cat_4xmmdr13s = 'IX/69/xmm4d13s'

    cat_chandra = 'IX/70/csc21mas'

    if DEBUG:
        help(Vizier)
        help(Vizier.get_catalog_metadata)
    # cat1 = cat_gaia_dr3
    cat1 = 'VII/294'
    catalog_metadata = Vizier(catalog=cat1).get_catalog_metadata()
    if DEBUG:
        help(catalog_metadata)

    catalog_metadata.info()
    #cat1 = cat_MQ
    cat1 = cat_chandra
    cat1 = cat_4xmmdr10st

    cat2 = cat_kids_dr5
    #Vizier(catalog=cat2).get_catalog_metadata()


    colname_xmatch_separation = 'angDist'
    outfile_xmatch = f"{table1_label}_xmatch_{table2_label}.fits"
    table_xmatch = cds_vizier_xmatch(
        table1=cat1, table2=cat2,
        max_distance=max_distance,
        colname_xmatch_separation=colname_xmatch_separation,
        outfile_xmatch=oufile_xmatch)
    print('Elapsed time:', time.time() - t0, 'seconds')
    table_xmatch.info()

    xdata = table_xmatch[colname_xmatch_separation]
    hist_binsize = 0.1
    bins = int(max_distance/hist_binsize)
    print(f'Number of bins, binsize: {bins} {hist_binsize}')
    plt.hist(xdata, bins=bins)
    plt.show()

    outfile = 'tmp_xmatch.fits'
    print(f'Write: {outfile}')
    table_xmatch.write(outfile, overwrite=True)
    print(f'Number of rows: {len(table_xmatch)}')
    print('Elapsed time:', time.time() - t0, 'seconds')

    logger.info('')

    return



if __name__ == "__main__":

    import fits_sql_history
    #help(fits_sql_history)

    logger = get_logger()
    logger.info('')
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Script: {sys.argv[0]}")
    logger.info(f"Arguments: {sys.argv[1:]}")

    # Connect to VizieR TAP service
    print(f'Use pyvo.dal.TAPService')
    tap_service = \
        pyvo.dal.TAPService("http://tapvizier.cds.unistra.fr/TAPVizieR/tap")


    tap_get_tile_list(tap_service_name=None,
                      table_name='II/383/kids_dr5',
                      colname='Tile')


    # test
    TEST = True
    if TEST:
        tap_get_radec_limits(tap_service_name=None,
                             tap_table='II/383/kids_dr5',
                             colnames_radec=['RAJ2000', 'DEJ2000'])

        sys.exit()


        logger = get_logger(test=False)

        tap_table_check(catalog="VII/289",
                        table="VII/289/dr16q",
                        table_colnames_radec=['RAJ2000', 'DEJ2000'],
                        radec_centre=[180.0, 0.0],
                        box_width_deg=0.1,
                        test=True,
                        verbose=True)

        sys.exit()

    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
    filename_timestamp = time.strftime('%Y%m%dT%H%M', time.gmtime())

    username = getpass.getuser()
    print('__name__:', __name__)

    logging.getLogger(__name__)

    logfile = 'kids_xmatch_cds_' + filename_timestamp + '.log'

    # set up file and screen logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s.%(msecs)03d %(name)-12s %(module)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
        datefmt='%y-%m-%dT%H:%M:%S',
        filename=logfile,
        filemode='w')

    # define a Handler which writes INFO messages or higher to sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d %(name)-12s: %(module)s %(funcName)s %(lineno)s %(levelname)-8s %(message)s',
        datefmt='%y-%m-%dT%H:%M:%S')
    # tell the handler to use this format
    console.setFormatter(formatter)

    # add the handler to the root logger
    logging.getLogger().addHandler(console)

    logging.info('Username: ' + username)
    logging.info(__name__)
    logging.info(__file__)

    logger = logging.getLogger(__name__)

    logger.info(f'Logfile created: {logfile}')

    t0 = time.time()

    main()

    sys.exit()

    from astroquery.utils.tap.core import TapPlus

    tap = TapPlus(url="https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
    job = tap.launch_job("SELECT column_name, unit, ucd "
                     "FROM TAP_SCHEMA.columns "
                     "WHERE table_name = 'VII/289/dr16q' "
                     "AND (ucd LIKE 'pos.eq.ra%' OR ucd LIKE 'pos.eq.dec%')")
    print(job.get_results())

    from astroquery.vizier import Vizier
    v = Vizier(columns=["**"])
    # Just fetch one row to inspect the table structure
    result = v.query_constraints(catalog="VII/289", RAJ2000=">0", DEJ2000=">0",
                              row_limit=1)
    print(result)
    for table in result:
        print()
        print(table.colnames)
        table.info()
        for icol, colname in enumerate(table.colnames):
            print(icol, colname)
            print(icol, table[colname].dtype)
            print(icol, table[colname].unit)
            print(icol, table[colname].info.description)
            if np.issubdtype(table[colname].dtype, np.number):
                print(f"Min: {np.min(table[colname])}")
                print(f"Max: {np.max(table[colname])}")
                print(f"Range: {np.max(table[colname]) - np.min(table[colname])}")
            else:
                # deal with masked array
                col = table[colname].data.data.astype(object)
                print(f"Min: {np.amin(col)}")
                print(f"Max: {np.amax(col)}")
            print()

    import pyvo
    import numpy as np

    tap = pyvo.dal.TAPService("https://tapvizier.cds.unistra.fr/TAPVizieR/tap")

    # Fetch a small sample - 3C 273 is RA=187.2779, Dec=+2.0524 (J2000)
    # It's one of the brightest quasars and definitely in DR16Q
    result = tap.search("""
        SELECT TOP 20 RAJ2000, DEJ2000, SDSS
        FROM "VII/289/dr16q"
        WHERE SDSS LIKE 'J1229%'
    """).to_table()

    print(result)

    # Check what scaling factor would recover sensible degrees
    for row in result:
        ra_raw = row['RAJ2000']
        dec_raw = row['DEJ2000']
        print(f"Raw: RA={ra_raw}, Dec={dec_raw}")
        for scale in [1, 1e3, 1e4, 3.6e6, 1e6, 1e7]:
            print(f"  /1e{np.log10(scale):.0f} -> RA={ra_raw/scale:.5f}, Dec={dec_raw/scale:.5f}")
     #3C 273 has SDSS name J122906.7+020308, so any row matching J1229 should have RA≈187.28°,
     # Dec≈+2.05°. The ratio of the raw value to those known coordinates gives you the exact
     # scale factor, and you can then hardcode it confidently in your cross-match query.
     # Sonnet 4.6


    sys.exit()
