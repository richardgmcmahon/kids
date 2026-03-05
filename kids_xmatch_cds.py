"""
https://cds-astro.github.io/tutorials/1_Intro_to_CDS_services_in_notebooks.html

https://astroquery.readthedocs.io/en/latest/api/astroquery.xmatch.XMatchClass.html


https://docs.g-vo.org/adql/

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

import fitsio
from fitsio import FITS,FITSHDR

DEBUG = False
if DEBUG:
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

def get_patch_limits_kids():

    dec_limits = [-40.0, 10.0]

    return


def tap_count_rows(table=None, columnName='RAJ2000'):
    """

    """
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

    return

def cds_vizier_xmatch(table1=None, table2=None,
                      colname_xmatch_separation = 'angDist',
                      max_distance=4.0,):

    cat1 = table1
    cat2 = table2

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
                          get_only_metadata=False,
                          xmatch_table_stats=True,
                          outfile_xmatch='cds_vizier_xmatch.fits',
                          outfile_xmatch_suffix=None,
                          add_nrows=True):
    """


    """

    import time
    from datetime import datetime, timezone

    import pyvo

    t0 = time.time()

    # Connect to CDS VizieR TAP service
    service = pyvo.dal.TAPService(
        "http://tapvizier.u-strasbg.fr/TAPVizieR/tap")
    service.describe()

    # print([tab_name for tab_name in service.tables.keys()])
    # help(service)

    # print(service.description)

    # Get table1 metadata
    if table1_name is None:
        table1_name = 'II/383/kids_dr5'

    print(f'table1_catalog: {table1_catalog}')

    GET_METADATA = False
    if GET_METADATA :
        print(f'Vizier(catalog=table1_catalog).get_catalog_metadata()')
        catalog_metadata = Vizier(catalog=table1_catalog).get_catalog_metadata()
        catalog_metadata.info()
        print(catalog_metadata)

    query = f"SELECT TOP 1 * FROM \"{table1_name}\""
    print(query)
    table1_info = service.search(query)

    print(f'table1_name: {table1_name}')
    print('table1_info.to_table().colnames')
    print(table1_info.to_table().colnames)

    if table2_name is None:
        table2_name = 'II/383/kids_dr5'
        table2_catalog = 'II/383'

    if GET_METADATA:
        print(f'Vizier(catalog=table2_catalog).get_catalog_metadata()')
        catalog_metadata = Vizier(catalog=table2_catalog).get_catalog_metadata()
        catalog_metadata.info()
        print(f'table2 catalog_metadata: {catalog_metadata}')

    query = f"SELECT TOP 1 * FROM \"{table2_name}\""
    print(query)
    table2_info = service.search(query)

    print(f'table2_name: {table2_name}')
    print('table2_info.to_table().colnames')
    print(table2_info.to_table().colnames)

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
          CIRCLE('ICRS', t2.{table2_colnames_radec[0]}, t2.{table2_colnames_radec[1]}, {max_distance}/3600.0)
     )
    ORDER BY separation_deg
    """

    print(f'ADQL query: {adql_query}')
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

    # help(service.search)
    print('Elapsed time:', time.time() - t0, 'seconds')
    results = service.search(adql_query)
    nrows = len(results)
    print(f'Query finished {nrows}')
    print('Elapsed time:', time.time() - t0, 'seconds')

    # Convert to astropy table
    table = results.to_table()
    print(f'Query results converted to table')
    print('Elapsed time:', time.time() - t0, 'seconds')

    # Optionally save to file as vot
    #table.write('milliquas_kids_crossmatch.vot', format='votable',
    #            overwrite=True)

    table.info()
    RUN_table_info_stats = False
    if RUN_table_info_stats:
        table.info(['attributes', 'stats'])

    print("\ncolumn dtypes:")
    print(table.dtype)


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

    print("\nUpdated dtypes:")
    print(table.dtype)

    table.info()


    if xmatch_table_stats:
        table.info(['attributes', 'stats'])

    # Display results
    print(f"\nFound {len(table)} matches within 10 arcseconds")
    table.pprint(max_lines=20, max_width=120)

    # fix the bad column u
    for col in table.colnames:
        try:
            # Test if this column's unit is FITS-compatible
            from astropy.io.fits import column
            if hasattr(table[col], 'unit') and table[col].unit is not None:
                str(table[col].unit.to_string('fits'))
        except Exception:
            print(f"Removing non-standard unit from column: {col}")
            table[col].unit = None


    # Optionally save to file
    print(f'outfile_xmatch template: {outfile_xmatch}')
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
    table.write(outfile, overwrite=True)
    print()
    print(f'{outfile} written successfully')
    print('Table write elapsed time:', time.time() - t1, 'seconds')
    print('Total elapsed time:', time.time() - t0, 'seconds')
    print()
    print(f'table1_name: {table1_name}')
    print(f'table2_name: {table2_name}')
    print(f'{len(table)} rows')
    print(f'max_distance (arcsec): {max_distance}')

    colname_xmatch_separation = 'separation_deg'
    itest = table[colname_xmatch_separation] < (max_distance/(2.0*3600.0))
    ntest = len(table[itest])
    print(len(table[itest]))

    itest = table[colname_xmatch_separation] < (max_distance/(4.0*3600.0))
    ntest = len(table[itest])
    print(len(table[itest]))



    # try fitsio option
    t1 = time.time()
    # fitsio.write(outputfile_xmatch, table.as_array(), clobber=True)
    # table.write(fitsio.write(outputfile_xmatch, clobber=True))
    print('Table write elapsed time:', time.time() - t1, 'seconds')
    print('Total elapsed time:', time.time() - t0, 'seconds')
    print()


    return table


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



def xmatch_input_tables(itable1=0,
                        itable2=1,
                        table1_label='csc12',
                        table2_label='kidsdr5'):

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


    table_name_list = ['VII/294/catalog',   # 0 MilliquasV8
                       'IX/70/csc21mas',    # 1 CSC21
                       'IX/69/xmm4d13s',    # 2 4XMMDR13
                       'IX/66/xmm411st',    # 3 4XMMDR11st
                       'J/ApJS/255/30/comp',# 4 VLASS_QL1
                       'J/A+A/661/A1/main', # 5
                       'J/A+A/661/A3/ctpmain',# 6
                       'VII/289/dr16q',     # 7
                       'I/356/qsocand',     # 8
                       'II/383/kids_dr5']   # 9 kids_dr5


    table_label_list = ['MQv8',
                        'csc21',
                        '4xmmdr13s',
                        '4xmmdr11st',
                        'VLASS_QL1',
                        '',
                        '',
                        'sdss_dr16q',
                        'gaiadr3_qsocand',
                        'kids_dr5']


    table_catalog_list =['VII/294',
                         'IX/70',
                         'IX/69',
                         'IX/66',
                         'J/ApJS/255/30',
                         'J/A+A/661/A1',
                         'J/A+A/661/A3',
                         'VII/289',
                         'I/356', #8
                         'II/383']   # 9


    table_colnames_radec = [['RAJ2000', 'DEJ2000'],
                            ['RAJ2000', 'DEJ2000']]

    table_name = ['name1', 'name2']
    table_catalog = ['catalog1', 'catalog2']
    table_label = [table_label_list[0],
    table_label_list[1]]



    table_name[1] = 'II/383/kids_dr5'
    table_catalog[1] = 'II/383'

    table_name[1] = 'VII/289/dr16q'
    table_catalog[1] = 'VII/289'

    table_catalog[0] = table_catalog_list[itable1]
    table_label[0] = table_label_list[itable1]
    table_name[0] = table_name_list[itable1]


    table_catalog[1] = table_catalog_list[itable2]
    table_label[1] = table_label_list[itable2]
    table_name[1] = table_name_list[itable2]



    if itable1 == 1:
        table_colnames_radec[0] = ['RAICRS', 'DEICRS']
    if itable2 == 1:
        table_colnames_radec[1] = ['RAICRS', 'DEICRS']

    if itable1 == 2:
        table_colnames_radec[0] = ['RA_ICRS', 'DE_ICRS']
    if itable2 == 2:
        table_colnames_radec[1] = ['RA_ICRS', 'DE_ICRS']

    if itable2 == 8:
        table_colnames_radec[1] = ['RA_ICRS', 'DE_ICRS']


    print(f'itable1: {itable1}')
    print(f'itable2: {itable2}')
    print(f'table_name: {table_name}')
    print(f'table_catalog: {table_catalog}')
    print(f'table_label: {table_label}')
    print(f'table_colnames_radec: {table_colnames_radec}')

    return table_name, table_catalog, table_label, table_colnames_radec


def main(DEBUG=False, tap_catalog=None, RUN_COUNT_ROWS=False):
    """


    """
    from datetime import datetime, timezone

    t0 = time.time()

    now = datetime.now(timezone.utc)
    # 2026-02-28 13:45:00.1 UTC
    print(now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-5] + " UTC")
    print(f'Elapsed time: {(time.time() - t0):.3f} seconds; {now.strftime("%Y-%m-%d %H:%M:%S.%f")} UTC')

    print(f'Elapsed time: {(time.time() - t0):.3f} seconds; {now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]} UTC')

    print(
    f'Elapsed time: {(time.time() - t0):.3f} seconds; '
    f'{now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]} UTC')

    # max_distance = 1.0
    # max_distance = 2.0
    # max_distance = 4.0
    max_distance = 8.0
    # max_distance = 16.0
    # max_distance = 32.0
    # max_distance = 60.0
    # max_distance = 120.0

    outfile_xmatch_suffix = f"an{int(max_distance)}as"

    itable1 = 1
    itable2 = 7

    table_names, table_catalogs, table_labels, tables_colnames_radec = \
    xmatch_input_tables(itable1=itable1,
    itable2=itable2,
    table1_label='csc12',
    table2_label='kidsdr5')


    table1_label = table_labels[0]
    table2_label = table_labels[1]

    outfile_xmatch = f"{table1_label}_xmatch_{table2_label}.fits"

    table1_name = table_names[0]
    table2_name = table_names[1]

    table1_catalog = table_catalogs[0]
    table2_catalog = table_catalogs[1]

    table1_colnames_radec = tables_colnames_radec[0]
    table2_colnames_radec = tables_colnames_radec[1]

    print(f'Running cds_vizier_xmatch_tap')
    cds_vizier_xmatch_tap(table1_name=table1_name,
                          table1_catalog=table1_catalog,
                          table2_name=table2_name,
                          table2_catalog=table2_catalog,
                          max_distance=max_distance,
                          get_only_metadata=True)


    print(f'RUN_COUNT_ROWS: {RUN_COUNT_ROWS}')
    if RUN_COUNT_ROWS:
        t1 = time.time()
        print('Elapsed time:', time.time() - t0, 'seconds')
        tap_count_rows(table=cat1)
        print('Elapsed time:', time.time() - t1, 'seconds')

        t1 = time.time()
        print('Elapsed time:', time.time() - t0, 'seconds')
        tap_count_rows(table=cat2)
        print('Elapsed time:', time.time() - t1, 'seconds')

        sys.exit()


    cds_vizier_xmatch_tap(table1_name=table1_name,
                          table1_catalog=table1_catalog,
                          table1_colnames_radec=table1_colnames_radec,
                          table2_name=table2_name,
                          table2_catalog=table2_catalog,
                          table2_colnames_radec=table2_colnames_radec,
                          max_distance=max_distance,
                          outfile_xmatch=outfile_xmatch,
                          outfile_xmatch_suffix=outfile_xmatch_suffix,
                          get_only_metadata=False)

    sys.exit()

    table_list = XMatch.get_available_tables()
    print(f'Number of tables: {len(table_list)}')
    #help(table_list)
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
    table_xmatch = cds_vizier_xmatch(
        table1=cat1, table2=cat2,
        max_distance=max_distance,
        colname_xmatch_separation=colname_xmatch_separation)
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


if __name__ == "__main__":
    RUN_COUNT_ROWS = False
    main(RUN_COUNT_ROWS=RUN_COUNT_ROWS)
