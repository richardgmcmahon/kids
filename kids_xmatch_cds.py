"""
https://cds-astro.github.io/tutorials/1_Intro_to_CDS_services_in_notebooks.html

https://astroquery.readthedocs.io/en/latest/api/astroquery.xmatch.XMatchClass.html

"""
# Standard Library
from pathlib import Path
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

# Astropy
from astropy import units as u
from astropy.table import Table

# Access astronomical databases
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch

import pyvo
from astroquery.utils.tap.core import TapPlus


help(XMatch)

"""
get_available_tables(self, *, cache=True)

query(self, cat1, cat2, max_distance, *,
colRA1=None, colDec1=None, colRA2=None, colDec2=None,
area='allsky', cache=True, get_query_payload=False, **kwargs)

query_async(self, cat1, cat2, max_distance, *,
colRA1=None, colDec1=None, colRA2=None, colDec2=None, area='allsky',
cache=True, get_query_payload=False, **kwargs)

"""

def tap_count_rows(table=None):
    import pyvo
    # Connect to VizieR TAP service
    tap = pyvo.dal.TAPService("http://tapvizier.cds.unistra.fr/TAPVizieR/tap")

    # Execute COUNT query
    query = f'SELECT COUNT(*) AS num_rows FROM "{table}"'

    result = tap.search(query)
    num_rows = result['num_rows'][0]
    print(f"Number of rows: {num_rows}")

    tap = TapPlus(url="http://tapvizier.cds.unistra.fr/TAPVizieR/tap")
    query = f'SELECT COUNT(*) FROM "{cat1}"'
    print(f'query: {query}')
    result = tap.launch_job(query).get_results()
    print(cat1, result)


def cds_vizier_xmatch(table1=None, table2=None,
                      max_distance=4.0):

    cat1 = table1
    cat2 = table2

    colname_xmatch_separation = 'angDist'
    max_distance = max_distance * u.arcsec
    print(max_distance)
    t0 = time.time()

    print(f'cat1: {cat1}')
    print(f'cat2: {cat2}')
    print(f'Maximum xmatch distance: {max_distance}')
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
    print()
    # 16  - 4 = 12

    print(f'Number of matchs > {max_distance/2}: {n_all - ntest_2} {(n_all - ntest_2)/12}')
    print(f'Number of matchs > {max_distance/4}: {n_all - ntest_4}')

    return table_xmatch


def cds_vizier_xmatch_tap(table1_name=None, table2_name=None,
                          max_distance=4.0):
    """


    """
    import pyvo
    # Connect to CDS VizieR TAP service
    service = pyvo.dal.TAPService(
        "http://tapvizier.u-strasbg.fr/TAPVizieR/tap")

    # Get table metadata
    if table1_name is None:
        table1_name = 'II/383/kids_dr5'
    query = f"SELECT TOP 1 * FROM \"{table1_name}\""
    table1_info = service.search(query)
    print('table1_info.to_table().colnames')
    print(table1_info.to_table().colnames)


    if table2_name is None:
        table2_name = 'II/383/kids_dr5'
    query = f"SELECT TOP 1 * FROM \"{table2_name}\""
    table2_info = service.search(query)
    print('table2_info.to_table().colnames')
    print(table2_info.to_table().colnames)


    # ADQL query to crossmatch Milliquas with KIDS DR5 within 10 arcsecond
    adql_query = """
    SELECT TOP 100
    mq.Name AS milliquas_name,
    mq.RAJ2000 AS mq_ra,
    mq.DEJ2000 AS mq_dec,
    mq.z AS redshift,
    kids.RAJ2000 AS kids_ra,
    kids.DEJ2000 AS kids_dec,
    kids.rmagGAAP AS kids_mag,
    DISTANCE(
        POINT('ICRS', mq.RAJ2000, mq.DEJ2000),
        POINT('ICRS', kids.RAJ2000, kids.DEJ2000)
    ) AS separation_deg
    FROM
    "VII/294/catalog" AS mq
    JOIN "II/383/kids_dr5" AS kids
     ON 1=CONTAINS(
          POINT('ICRS', mq.RAJ2000, mq.DEJ2000),
          CIRCLE('ICRS', kids.RAJ2000, kids.DEJ2000, 10.0/3600.0)
     )
    ORDER BY separation_deg
    """

    print(f'ADQL query: {adql_query}')
    # Execute the query
    print("Running crossmatch query...")
    results = service.search(adql_query)

    # Convert to astropy table
    table = results.to_table()

    table.info()

    # Display results
    print(f"\nFound {len(table)} matches within 10 arcseconds")
    table.pprint(max_lines=20, max_width=120)

    # Optionally save to file
    table.write('milliquas_kids_crossmatch.fits', overwrite=True)


    sys.exit()


def fix_votable_object(table, dtype_new='bool',
                       verbose=False,
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



def main():


    t0 = time.time()

    table1_name = 'VII/294/catalog'
    table2_name = 'II/383/kids_dr5'
    cds_vizier_xmatch_tap(table1_name=table1_name,
                          table2_name=table2_name)

    table_list = XMatch.get_available_tables()
    print(f'Number of tables: {len(table_list)}')
    help(table_list)
    #print(table_list)

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

    help(Vizier)
    help(Vizier.get_catalog_metadata)
    # cat1 = cat_gaia_dr3
    cat1 = 'VII/294'
    catalog_metadata = Vizier(catalog=cat1).get_catalog_metadata()
    help(catalog_metadata)
    catalog_metadata.info()
    cat1 = cat_MQ
    cat2 = cat_kids_dr5
    #Vizier(catalog=cat2).get_catalog_metadata()


    RUN_COUNT_ROWS = False
    if RUN_COUNT_ROWS:
        t1 = time.time()
        print('Elapsed time:', time.time() - t0, 'seconds')
        tap_count_rows(table=cat1)
        print('Elapsed time:', time.time() - t1, 'seconds')

        t1 = time.time()
        print('Elapsed time:', time.time() - t0, 'seconds')
        tap_count_rows(table=cat2)
        print('Elapsed time:', time.time() - t1, 'seconds')

    max_distance = 8.0
    max_distance = 16.0

    table_xmatch = cds_vizier_xmatch(
        table1=cat1, table2=cat2, max_distance=max_distance)

    xdata = table_xmatch[colname_xmatch_separation]
    hist_binsize = 0.1
    bins = max_distance/hist_binsize
    print(f'Number of bins, binsize: {bins} {hist_binsize}')
    plt.hist(xdata, bins=bins)
    plt.show()

    outfile = 'tmp_xmatch.fits'
    print(f'Write: {outfile}')
    table_xmatch.write(outfile, overwrite=True)
    print(f'Number of rows: {len(table_xmatch)}')
    print('Elapsed time:', time.time() - t0, 'seconds')


if __name__ == "__main__":
    main()
