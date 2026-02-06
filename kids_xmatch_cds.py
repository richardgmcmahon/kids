"""
https://cds-astro.github.io/tutorials/1_Intro_to_CDS_services_in_notebooks.html

https://astroquery.readthedocs.io/en/latest/api/astroquery.xmatch.XMatchClass.html

"""
# Standard Library
from pathlib import Path
import os
import time

# Astropy
from astropy import units as u
from astropy.table import Table

# Access astronomical databases
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch



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
print(Vizier(catalog=cat1).get_catalog_metadata())
cat1 = cat_MQ

cat2 = cat_kids_dr5
#Vizier(catalog=cat2).get_catalog_metadata()

from astroquery.utils.tap.core import TapPlus

tap = TapPlus(url="http://tapvizier.cds.unistra.fr/TAPVizieR/tap")
query = "SELECT COUNT(*) FROM \"VII/294/catalog\""
print(f'query: {query}')
query = f'SELECT COUNT(*) FROM "{cat1}"'
result = tap.launch_job(query).get_results()
print(cat1, result)

query = f'SELECT COUNT(*) FROM "{cat2}"'
result = tap.launch_job(query).get_results()
print(cat2, result)


max_distance = 5.0
max_distance = max_distance * u.arcsec
t0 = time.time()
table_xmatch = XMatch.query(cat1, cat2, max_distance)
print('Elapsed time:', time.time() - t0, 'seconds')

print(f'Number of rows: {len(table_xmatch)}')
table_xmatch.info()
table_xmatch.info(['attributes', 'stats'])

outfile = 'tmp_xmatch.fits'
print('Write: {outfile}')
table_xmatch.write(outfile, overwrite=True)
print(f'Number of rows: {len(table_xmatch)}')
print('Elapsed time:', time.time() - t0, 'seconds')
