"""

https://docs.astropy.org/en/stable/api/astropy.table.JSViewer.html#astropy.table.JSViewer

Lai+2024

Sample of 83 quasars with 4.5 < z < 5.3

https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.3912L
https://arxiv.org/abs/2311.05098

https://github.com/samlaihei/XQz5

See also:
https://ui.adsabs.harvard.edu/abs/2016ApJ...819...24W
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/819/24
https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/ApJ/819/24/table1

Datasets

LS DR10

https://www.legacysurvey.org/dr10/description/#photometry


SDSS DR9
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=V/139

http://tapvizier.cds.unistra.fr/adql/?

PS1:

unWISE: https://ui.adsabs.harvard.edu/abs/2019ApJS..240...30S
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/363

CATWISE
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/365

MilliQuas
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=VII/290
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=VII/294

Wang+2016
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/819/24

144MHz radio fluxes of z>5 quasars (Gloudemans+, 2021)
2021A&A...656A.137G
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A%2BA/656/A137

Unveiling weak radio quasar population at z≥4 (Perger+, 2019)
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/490/2542
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/other/FrASS/4.9

Ross&Cross+2020
488 very high-z quasars with near and MIR photometry (
https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/494/789

https://datalab.noirlab.edu/xmatch.php

PS1 data TAP interface
https://spacetelescope.github.io/notebooks/notebooks/MAST/PanSTARRS/PS1_DR2_TAP/PS1_DR2_TAP.html

Fast generation of FITS images cutouts from HiPS datasets
https://alasky.cds.unistra.fr/hips-image-services/hips2fits
https://astroquery.readthedocs.io/en/latest/hips2fits/hips2fits.html


https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/convcoord/convcoord.pl



"""
import os
import sys
import time

import numpy as np
import numpy.ma as ma
from scipy import stats

import matplotlib.pyplot as plt

import pandas as pd

import astropy
print('astropy.__version__:', astropy.__version__)
#from astropy.metadata import version
import astropy.version
#help(astropy.version)
print('astropy.__version__:', astropy.version)

from astropy.table import Column, Table, hstack
from astropy import units as u
from astropy.coordinates import SkyCoord

import xmatch as xm
#from xmatch_ls import *

#help(xm)
#help(xm.xmatch_cat)
#help(xm.xmatch_checkplots)

from lineInfo import *
import rgm_util.mk_urls as mkurl
#import mk_urls as mkurl
help(mkurl)

t0 = time.time()


def LS_Flux_to_Mag(flux=None):
    """https://www.legacysurvey.org/dr10/description/#photometry

    The brightnesses of objects are all stored as linear fluxes in units of
    nanomaggies. The conversion from linear fluxes to magnitudes is
    𝑚=22.5−2.5log10(flux). These linear fluxes are well-defined even at the
    faint end, and the errors on the linear fluxes should be very close to a
    normal distribution. The fluxes can be negative for faint objects, and
    indeed we expect many such cases for the faintest objects.

    Linear flux units, where an object with an AB magnitude of 0 has a flux
    of 1.0 maggie. A convenient unit is the nanomaggie: a flux of 1
    nanomaggie corresponds to an AB magnitude of 22.5.

    nanomaggie

    Linear flux units, where an object with an AB magnitude of 22.5 has a
    flux of 1×10−9 maggie or 1.0 nanomaggie.

    10^-9 = 22.5

    """

    # NaNs are returned where x is negative.
    mag = 22.5 - (2.50 * np.log10(flux))

    return mag


def invar_to_sigma(invar=None):

    sigma = np.sqrt(1.0/invar)

    return sigma




def read_xqz5(path='./', filename=None, fmt=None, fix=False):

    infile = path + filename

    print('Reading:', infile)
    table = Table.read(infile)
    print(len(table), 'rows read in')
    table.info()
    table.info(['attributes', 'stats'])

    # ID is object name
    print(table['ID'])
    # PID is ESO Programme ID
    print(table['PID'])

    # replace msun with Msun which is a FITs standard unit
    # help(table['logMbh'])
    #print(table['logMbh'].unit)
    table['logMbh'].unit = 'Msun'
    table['elogMbh'].unit = 'Msun'
    #print(table['logMbh'].unit)

    #sys.exit()

    # Source is Telescope data source
    print(table['Source'])

    # save RA, Dec in csv file for cross matching

    WRITE_CSV = False
    if WRITE_CSV:
        colnames_radec = ('RA', 'Dec')

        out_table = table['RA', 'Dec']
        out_table.info()
        out_table.write('tmp.csv', overwrite=True)
        out_table.write('tmp.ecsv', overwrite=True)

        # write out whole file as csv
        outfilename = os.path.splitext(filename)[0] + '.csv'
        print('Writing:', outfilename)
        table.write(outfilename, overwrite=False)

    fix_RADEC= False
    if fix_RADEC is True:
        fix_name = 'SDSS J095727.86+051905.2'
        print('Fix the RA, Dec for', fix_name)
        # using DR16Q Lyke et al., 2020
        RA_fixed = 149.366112
        Dec_fixed = 5.318132

        fix_name = 'SDSS J095712.20+101618.5'
        # using https://ui.adsabs.harvard.edu/abs/2017AJ....153..184Y
        print('Fix the RA, Dec for', fix_name)
        RA_fixed = 149.300833
        Dec_fixed = 10.271806

    # Check the positions against the RADEC names:
    check_radec = True
    if check_radec:

        nrows = len(table)
        ra  = np.empty([nrows], dtype=float)
        dec  = np.empty([nrows], dtype=float)
        sep_arcsec  = np.empty([nrows], dtype=float)
        for i, row in enumerate(table):
            radec_string = table[i]['ID']
            ra[i], dec[i] =  mkurl.IAUname_to_radec(
                radec_string=radec_string)
            c1 = SkyCoord(table[i]['RA']*u.degree, table[i]['Dec']*u.degree,
                          frame='icrs')
            c2 = SkyCoord(ra[i]*u.degree, dec[i]*u.degree,
                          frame='icrs')
            sep_arcsec[i] = c1.separation(c2).arcsec

            print(i, table[i]['ID'],
                  ra[i], dec[i],
                  table[i]['RA'], table[i]['Dec'],
                  (table[i]['RA'] - ra[i]) * 3600.0,
                  (table[i]['Dec'] - dec[i]) * 3600.0,
                  sep_arcsec[i])

        print()
        for i, row in enumerate(table):
            print(i, table[i]['ID'],
                  table[i]['RA'], table[i]['Dec'],
                  sep_arcsec[i])

    return table



def read_data(path='./', filename=None, fmt=None):

    infile = path + filename

    print('Reading:', infile)
    table = Table.read(infile)
    print(len(table), 'rows read in')
    table.info()
    table.info(['attributes', 'stats'])

    if fmt == 'xqz5':
        # ID is object name
        #print(table['ID'])
        # PID is ESO Programme ID
        #print(table['PID'])

        # Source is Telescope
        print(table['Source'])

    # save RA, Dec in csv file for cross matching

    WRITE_CSV = False
    if WRITE_CSV:
        colnames_radec = ('RA', 'Dec')

        out_table = table['RA', 'Dec']
        out_table.info()
        out_table.write('tmp.csv', overwrite=True)
        out_table.write('tmp.ecsv', overwrite=True)


    return table



def plot_redshift_mag(table=None, filename=None):
    """

    """


    # plot mag versus redshift locus
    plt.figure(figsize=(8.0, 6.0))

    xdata = table['MgII Redshift']
    ydata = table['i_mag']

    ndata = len(xdata)
    label = str(ndata)
    plt.title(filename)
    plt.scatter(xdata, ydata, label=label)

    plt.xlim(4.40, 5.40)
    plt.ylim(16.5, 21.5)

    plt.gca().invert_yaxis()
    plt.legend(loc='upper right')

    plt.xlabel('Redshift (MgII)')
    plt.ylabel('i mag')


    plt.show()

    #
    plt.figure(figsize=(8.0, 6.0))
    ydata = table['z_mag']

    ndata = len(xdata)
    label = str(ndata)
    plt.title(filename)
    plt.scatter(xdata, ydata, label=label)

    plt.xlim(4.40, 5.40)
    plt.ylim(16.5, 21.5)

    plt.gca().invert_yaxis()
    plt.legend(loc='upper right')

    plt.ylabel('z mag')

    plt.show()

    #
    plt.figure(figsize=(8.0, 6.0))
    ydata = table['Hmag']

    ndata = len(xdata)
    label = str(ndata)
    plt.title(filename)
    plt.scatter(xdata, ydata, label=label)

    plt.xlim(4.40, 5.40)
    plt.ylim(14.5, 19.5)

    plt.gca().invert_yaxis()
    plt.legend(loc='upper right')

    plt.xlabel('Redshift (MgII)')
    plt.ylabel('H mag')

    plt.show()

    return


def explore_ls(table=None):

    table.info(['attributes', 'stats'])
    print(len(table), 'rows')

    return


def explore_PS1(table=None, stack=True,
                  verbose=False, debug=False):
    """
    when there is no match all the right hand columns are 0.0

    """

    table.info(['attributes', 'stats'])
    print(len(table), 'rows')

    if stack:
        itest = (table['primaryDetection'] == 0)
        print('primaryDetection = 0:',
          len(table[itest]), len(table))

        itest = (table['primaryDetection'] == 1)
        print('primaryDetection = 1:',
          len(table[itest]), len(table))


        itest = (table['primaryDetection'] > 1)
        print('primaryDetection > 1:',
          len(table[itest]), len(table))

        itest = (table['primaryDetection'] < 0)
        print('primaryDetection < 1:',
          len(table[itest]), len(table))


    itest = (table['_searchID_'] > 0)
    print('Number of valid  _searchID_:',
          len(table[itest]), len(table))

    unique =  np.unique(table['_searchID_'])
    print('Number of unique _searchID_ :', len(unique), len(table))

    unique =  np.unique(table['_ra_'])
    print('Number of unique _ra_ :', len(unique), len(table))

    unique =  np.unique(table['_dec_'])
    print('Number of unique _dec_ :', len(unique), len(table))

    unique =  np.unique(table['MatchRA'])
    print('Number of unique MatchRA :', len(unique), len(table))

    unique =  np.unique(table['MatchDEC'])
    print('Number of unique MatchDEC :', len(unique), len(table))

    unique =  np.unique(table['MatchID'])
    print('Number of unique MatchID:', len(unique), len(table))

    unique =  np.unique(table['objName'])
    print('Number of unique objName:', len(unique), len(table))

    isort = np.argsort(table['_searchID_'])

    for i, row in enumerate(table[isort]):
        print(i, row['_searchID_'], row['MatchID'],
              row['dstArcSec'],
              row['objName'],
              '{:.5f}'.format(row['_ra_']),
              '{:.5f}'.format(row['_dec_']),
              '{:.5f}'.format(row['MatchRA']),
              '{:.5f}'.format(row['MatchDEC']),
              row['primaryDetection'])



    return



def explore_UHS(table=None, verbose=False, debug=False):
    """
    when there is no match all the right hand columns are 0.0

    """

    table.info(['attributes', 'stats'])
    print(len(table), 'rows')

    itest = (table['SOURCEID'] > 0)
    print('Number of valid  SOURCEID rows:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) | (table['KAPERMAG3'] > 1.0)
    print('Number of valid rows:', len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0)
    print('Number of valid  JAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['KAPERMAG3'] > 1.0)
    print('Number of valid  KAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) & (table['KAPERMAG3'] > 1.0)
    print('Number of sources delected in J AND K bands:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) | (table['KAPERMAG3'] > 1.0)
    print('Number of sources delected in J OR K bands:',
          len(table[itest]), len(table))

    if verbose and debug:
        for irow, row in enumerate(table):
            print(irow, row)

    return



def explore_UKIDSS_LAS(table=None, verbose=False, debug=False):
    """
    when there is no match all the right hand columns are 0.0

    """

    table.info(['attributes', 'stats'])
    print(len(table), 'rows')

    itest = (table['SOURCEID'] > 0)
    print('Number of valid  SOURCEID rows:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) | (table['KAPERMAG3'] > 1.0)
    print('Number of valid rows:', len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0)
    print('Number of valid  JAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['KAPERMAG3'] > 1.0)
    print('Number of valid  KAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) & (table['KAPERMAG3'] > 1.0)
    print('Number of sources delected in J AND K bands:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) | (table['KAPERMAG3'] > 1.0)
    print('Number of sources delected in J OR K bands:',
          len(table[itest]), len(table))

    if verbose and debug:
        for irow, row in enumerate(table):
            print(irow, row)

    return



def explore_VISTA(table=None, vhs=True, viking=False,
                  verbose=False, debug=False):
    """
    when there is no match all the right hand columns are 0.0

    """

    table.info(['attributes', 'stats'])
    print(len(table), 'rows')

    itest = (table['SOURCEID'] > 0)
    print('Number of valid  SOURCEID rows:',
          len(table[itest]), len(table))

    if vhs:
        itest = (table['YAPERMAG3'] > 1.0) | \
            (table['JAPERMAG3'] > 1.0) | \
            (table['HAPERMAG3'] > 1.0) | \
            (table['KSAPERMAG3'] > 1.0)
        print('Number of valid rows:', len(table[itest]), len(table))

    if viking:
        itest = \
            (table['ZAPERMAG3'] > 1.0) | \
            (table['YAPERMAG3'] > 1.0) | \
            (table['JAPERMAG3'] > 1.0) | \
            (table['HAPERMAG3'] > 1.0) | \
            (table['KSAPERMAG3'] > 1.0)
        print('Number of valid rows:', len(table[itest]), len(table))

    if viking:
        itest = (table['ZAPERMAG3'] > 1.0)
        print('Number of valid  ZAPERMAG3 rows:',
              len(table[itest]), len(table))


    itest = (table['YAPERMAG3'] > 1.0)
    print('Number of valid  YAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0)
    print('Number of valid  JAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['HAPERMAG3'] > 1.0)
    print('Number of valid  HAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['KSAPERMAG3'] > 1.0)
    print('Number of valid  KSAPERMAG3 rows:',
          len(table[itest]), len(table))

    itest = (table['YAPERMAG3'] > 1.0) & (table['JAPERMAG3'] > 1.0) & \
        (table['HAPERMAG3'] > 1.0) & (table['KSAPERMAG3'] > 1.0)
    print('Number of sources delected in all bands:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) & (table['KSAPERMAG3'] > 1.0)
    print('Number of sources delected in J AND Ks bands:',
          len(table[itest]), len(table))

    itest = (table['JAPERMAG3'] > 1.0) | (table['KSAPERMAG3'] > 1.0)
    print('Number of sources delected in J OR Ks bands:',
          len(table[itest]), len(table))

    if verbose and debug:
        for irow, row in enumerate(table):
            print(irow, row)

    return



def xmatch_ls(table_xqz5=None, table_ls=None,
              title_prefix='None',
              suptitle_prefix='None',
              plotfile_prefix=None,
              checkplots=True):
    """



    """

    ra1 = table_xqz5['RA']
    dec1 = table_xqz5['Dec']

    ra2 = table_ls['ra']
    dec2 = table_ls['dec']

    colnames_radec1 = ['RA', 'Dec']
    colnames_radec2 = ['ra', 'dec']

    table1 = table_xqz5
    table2 = table_ls
    # returns (idx1, idx2), dr, dra, ddec
    idx2, dr, dra, ddec = xm.xmatch_cat(table1=table1,
                                        table2=table2,
                                        colnames_radec1=colnames_radec1,
                                        colnames_radec2=colnames_radec2,
                                        verbose=True)
    print('Elapsed time(secs):', time.time() - t0)
    print('dr range:', np.min(dr), np.max(dr))
    print('dra range:', np.min(dra), np.max(dra))
    print('dddec range:', np.min(ddec), np.max(ddec))
    print()

    itest = (dr < 1.5)
    itest_bad = (dr > 1.5)
    print('itest_bad xmatch:', itest_bad)
    idx2 = idx2[itest]
    dr = dr[itest]
    table_xqz5 = table_xqz5[itest]
    table_ls = table_ls[idx2]
    print()

    table1 = table_xqz5
    table2 = table_ls
    for i, idx in enumerate(idx2):
       print(i, idx,
             dr[i], table1['ID'][i],
             table1[colnames_radec1[0]][i],
             table1[colnames_radec1[1]][i],
             table2[colnames_radec2[0]][idx],
             table2[colnames_radec2[1]][idx])

    if checkplots:
        print('Run xmatch_checkplots')
        ra1 = table_xqz5['RA']
        dec1 = table_xqz5['Dec']
        ra2 = table_ls['ra']
        dec2 = table_ls['dec']
        xm.xmatch_checkplots(ra1=ra1, dec1=dec1,
                             ra2=ra2, dec2=dec2,
                             rmax=5.0,
                             suptitle='XQz5 xm LS_DR10',
                             plotfile_prefix='XQz5_xm_LS_DR10_')
        print('Elapsed time(secs):', time.time() - t0)


    table_xmatch = hstack([table_xqz5, table_ls[idx2]])
    table_xmatch.info(['attributes', 'stats'])

    table_xmatch.write('xmatch.fits', overwrite=True)

    return table_xmatch

def mk_webpage(table=None,
               colname_name=None,
               colname_ra='RA',
               colname_dec='Dec',
               colname_ra_marker=None,
               colname_dec_marker=None,
               colname_list=None,
               gaia_edr3=True,
               colname_list_format=None):
    """

    # https://docs.astropy.org/en/stable/api/astropy.table.JSViewer.html#astropy.table.JSViewer
    # https://docs.astropy.org/en/stable/api/astropy.io.ascii.HTML.html#astropy.io.ascii.HTML

    # https://docs.astropy.org/en/stable/api/astropy.table.JSViewer.html

    table.show_in_notebook()

    table.show_in_browser()

    For a 5-sigma point source detection limit in g,
    5/(sqrt psfdepth_g) gives flux in nanomaggies and
    -2.5[log10(5 / (sqrt psfdepth_g)) - 9]
    gives corresponding AB magnitude

    psfdepth_mag_g =

    """

    nrows = len(table)

    # create empty table for html webpage
    table_out = Table()

    # Index column
    table_out['Id'] = np.linspace(1, nrows, nrows, dtype=int)

    if colname_name is None:
        table_out['Name'] = table['ID']
    if colname_name is not None:
        table_out['Name'] = table[colname_name]

    if colname_ra is None:
        table_out['RA'] = table['RA']
    if colname_ra is not None:
        table_out['RA'] = table[colname_ra]
    table_out['RA'].format = '{:10.5f}'

    if colname_dec is None:
        table_out['Dec'] = table['Dec']
    if colname_dec is not None:
        table_out['Dec'] = table[colname_dec]
    table_out['Dec'].format = '{:10.5f}'

    if colname_ra_marker is not None:
        RA_marker = table[colname_ra_marker]
    if colname_dec_marker is not None:
        Dec_marker = table[colname_dec_marker]

    if colname_list is not None:
        for icol, colname in enumerate(colname_list):
            table_out[colname] = table[colname]
            if colname_list_format is not None:
                if colname_list_format[icol] is not None:
                    table_out[icol].format = colname_list_format[icol]

    MK_THIN_TABLE = True
    if not MK_THIN_TABLE:
        table_out['dist'] = table['dist_arcsec']
        table_out['dist'].format = '{:6.3f}'


        table_out['MgII Redshift'] = table['MgII Redshift']
        table_out['MgII Redshift'].format = '{:8.4f}'

        table_out['eMgII Redshift'] = table['eMgII Redshift']
        table_out['eMgII Redshift'].format = '{:8.4f}'

        table_out['gaia_g'] = table['gaia_phot_g_mean_mag']
        table_out['gaia_g'].format = '{:6.2f}'

        table_out['type'] = table['type']

        table_out['mag_g'] = table['mag_g']
        table_out['mag_g'].format = '{:6.2f}'
        table_out['flux_g'] = table['flux_g']
        table_out['ivar_g'] = table['flux_ivar_g']
        table_out['sigma flux_g'] = invar_to_sigma(table['flux_ivar_g'])
        table_out['sigma mag_g'] = LS_Flux_to_Mag(table_out['sigma flux_g'])
        table_out['sigma mag_g'].format = '{:6.2f}'

        table_out['snr_g'] = table['snr_g']
        table_out['snr_g'].format = '{:6.2f}'
        table_out['psfdepth_g'] = table['psfdepth_g']

        # -2.5 [ log10 (5 / (sqrt psfdepth_g)) - 9]
        table_out['psfdepth mag_g'] = \
            -2.5 * (np.log10(5 / (np.sqrt(table['psfdepth_g']))) - 9)
        table_out['psfdepth mag_g'].format = '{:6.2f}'

        table_out['mag_r'] = table['mag_r']
        table_out['mag_r'].format = '{:6.2f}'
        table_out['flux_r'] = table['flux_r']
        table_out['ivar_r'] = table['flux_ivar_r']
        table_out['snr_r'] = table['snr_r']
        table_out['psfdepth_r'] = table['psfdepth_r']

        table_out['mag_i'] = table['mag_i']
        table_out['mag_i'].format = '{:6.2f}'
        table_out['flux_i'] = table['flux_i']
        table_out['ivar_i'] = table['flux_ivar_i']
        table_out['snr_i'] = table['snr_i']
        table_out['psfdepth_i'] = table['psfdepth_i']

        table_out['mag_z'] = table['mag_z']
        table_out['mag_z'].format = '{:6.2f}'
        table_out['flux_z'] = table['flux_z']
        table_out['ivar_z'] = table['flux_ivar_z']
        table_out['snr_z'] = table['snr_z']
        table_out['psfdepth_z'] = table['psfdepth_z']

        table_out['g_r'] = table['g_r']
        table_out['r_i'] = table['r_i']
        table_out['r_z'] = table['r_z']
        table_out['i_z'] = table['i_z']


    table_out.show_in_browser(
        jsviewer=True,
        table_class='compact cell-border hover')


    links_CDS_Simbad = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_CDS_Simbad(RA=RA, Dec=Dec, radius=3.0)
        # print(url_LS)
        links_CDS_Simbad[i] = '<a href="' + url + '"> Simbad </a>'

    links_CDS_Vizier = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_CDS_Vizier(RA=RA, Dec=Dec, radius=15.0)
        # print(url_LS)
        links_CDS_Vizier[i] = '<a href="' + url + '"> Vizier </a>'

    links_CDS_Vizier_GaiaDR3 = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_CDS_Vizier(catalogue='I/355/gaiadr3',
                                       RA=RA, Dec=Dec, radius=15.0)
        # print(url_LS)
        links_CDS_Vizier_GaiaDR3[i] = '<a href="' + url + '"> Gaia </a>'

    # KIDS DR5

    cds_vizier_table = 'II/383/kids_dr5'
    links_CDS_Vizier_KIDS_DR5 = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_CDS_Vizier(catalogue='II/383/kids_dr5',
                                       RA=RA, Dec=Dec, radius=15.0)
        # print(url_LS)
        links_CDS_Vizier_KIDS_DR5[i] = '<a href="' + url + '"> KIDS_DR5 </a>'


# Milliquas v8
    links_CDS_Vizier_MQv8 = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_CDS_Vizier(catalogue='VII/294',
                                       RA=RA, Dec=Dec, radius=15.0)
        # print(url_LS)
        links_CDS_Vizier_MQv8[i] = '<a href="' + url + '"> MQv8 </a>'


    LS_links = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        RA_marker = table[i][colname_ra_marker]
        Dec_marker = table[i][colname_dec_marker]
        print(i, RA, Dec)
        print(i, RA_marker, Dec_marker)
        url = mkurl.mk_url_LS(RA=RA, Dec=Dec, DR='dr10',
                              gaia_edr3=gaia_edr3,
                              RA_marker=RA_marker,
                              Dec_marker=Dec_marker)
        # print(url_LS)
        LS_links[i] = '<a href="' + url + '"> LS </a>'


    eROSITA_links = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_eROSITA(RA=RA, Dec=Dec)
        # print(url_LS)
        eROSITA_links[i] = '<a href="' + url + '"> eROSITA </a>'


    NED_links = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_NED(RA=RA, Dec=Dec, radius=10.0)
        # print(url_LS)
        NED_links[i] = '<a href="' + url + '"> NED </a>'


    SDSS_links = np.empty(nrows, dtype=object)
    for i, row in enumerate(table):
        RA = table[i][colname_ra]
        Dec = table[i][colname_dec]
        url = mkurl.mk_url_sdss_skyserver(
            RA=RA, Dec=Dec, DR='dr18')
        SDSS_links[i] = '<a href="' + url + '"> SDSS </a>'

    # use some blanks in the names so that web page column pages have
    # line breaks in column titles
    table_out.add_column(LS_links, index=4, name='LS DR10')
    table_out.add_column(SDSS_links, index=5, name='SDSS')
    table_out.add_column(NED_links, index=6, name='NED')
    table_out.add_column(eROSITA_links, index=6, name='eROSITA')
    table_out.add_column(links_CDS_Simbad, index=7,
                         name='CDS Simbad')
    table_out.add_column(links_CDS_Vizier, index=8,
                         name='CDS Vizier')
    table_out.add_column(links_CDS_Vizier_GaiaDR3, index=9,
                         name='Gaia DR3')
    table_out.add_column(links_CDS_Vizier_MQv8, index=10,
                         name='MQ v8')
    table_out.add_column(links_CDS_Vizier_KIDS_DR5, index=11,
                             name='KIDS_DR5')


    print('Write:', 'tmp_jsviewer.html')
    # blanks mean that webpage text wraps to data column widths
    htmldict = {'raw_html_cols': ['LS DR10', 'SDSS', 'NED', 'eROSITA',
                                  'CDS Simbad', 'CDS Vizier',
                                  'Gaia DR3', 'MQ v8', 'KIDS_DR5']}
    # uses https://datatables.net/
    # https://www.datatables.net/manual/styling/classes
    table_out.write('tmp_jsviewer.html',
                    format='jsviewer',
                    max_lines = 5000,
                    htmldict=htmldict,
                    # jskwargs=jskwargs,
                    table_class='compact cell-border hover',
                    overwrite=True)

    #LineInfo()

    return



if __name__=='__main__':


    # https://developer.lsst.io/stack/logging.html#basic-usage-in-python
    import logging

    logging.basicConfig(
        level=logging.INFO,
        format="{name} {asctime} {levelname}: {message}",
        style="{")

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    #logging.info('Some information during normal operation')
    #logging.warning('Here is a warning!')

    logger.info('Some information during normal operation')
    logger.warning('Here is a warning!')

    logger.info("This is information about the star selector algorithm execution. %f", 3.14)



    xmatch_checkplots = True

    READ_PS1 = False
    READ_UHS = False
    READ_UKIDSS_GCS = False
    READ_UKIDSS_LAS = False

    READ_VHS = False
    READ_VIKING = False

    PLOT_REDSHIFT = False

    READ_XQZ5 = True

    READ_XMM_XMATCH= True

    if READ_XMM_XMATCH:
        filename = 'test_table.fits'
        filename = 'worst_LR_most_likely_neighbor.fits'
        path = './'
        infile = path + filename
        print('Reading: {infile}')
        table = Table.read(infile)
        plotfile_prefix = filename
        table.info(['attributes', 'stats'])

        # create html webpage
        # colname_name='4XMM',
        colname_ra='RA_ICRS'
        colname_dec='DEC_ICRS'

        colname_name = 'source_id'
        colname_ra='RA_src'
        colname_dec='DEC_src'
        colname_ra_marker = 'RA_ctp'
        colname_dec_marker = 'DEC_ctp'

        colname_list=['mag', 'lr', 'd_ell', 'R_1SIG_MAJ', 'total_matches']
        mk_webpage(table=table,
                   colname_name=colname_name,
                   colname_ra=colname_ra,
                   colname_dec=colname_dec,
                   colname_ra_marker=colname_ra_marker,
                   colname_dec_marker=colname_dec_marker,
                   colname_list=colname_list)

        sys.exit()

    if READ_XQZ5:
        path = './'
        filename = 'XQz5_properties_20231228.fits'
        table_xqz5 = read_xqz5(filename=filename, path=path)
        table_xqz5.info(['attributes', 'stats'])
        plotfile_prefix = filename

        #sys.exit()

    if PLOT_REDSHIFT:
        plot_redshift_mag(table=table_xqz5, filename=filename)

        sys.exit()


    if READ_PS1:
        path = './XMatch/'
        filename = 'xqz5_xm_PS1DR2_Stack_an3as.fits'
        #filename = 'xqz5_xm_PS1DR2_Mean_an3as.fits'
        #filename = 'xqz5_xm_PS1DR2_ForcedMean_an3as.fits'
        infile = path + filename
        print('Reading:', infile)
        table_ps1 = Table.read(infile)
        print(len(table_ps1), 'rows')

        explore_PS1(table=table_ps1, stack=True)

        sys.exit()


    if READ_UHS:
        path_uhs = './XMatch/'
        filename_uhs = 'xqz5_xm_UHSDR2_nn5as.fits'
        infile_uhs = path_uhs + filename_uhs
        table_uhs = Table.read(infile_uhs)
        print(len(table_uhs), 'rows')

        explore_UHS(table=table_uhs)


    if READ_UKIDSS_LAS:
        path = './XMatch/'
        filename_ukidss = 'xqz5_xm_UKIDSSDR11PLUS_LAS_nn5as.fits'
        infile = path + filename_ukidss
        print('Reading:', infile)
        table_ukidss = Table.read(infile)
        print(len(table_ukidss), 'rows')

        explore_UKIDSS_LAS(table=table_ukidss)


    if READ_UKIDSS_GCS:
        path = './XMatch/'
        filename_ukidss = 'xqz5_xm_UKIDSSDR11PLUS_GCS_nn5as.fits'
        infile_ukidss = path + filename_ukidss
        table_ukidss = Table.read(infile_ukidss)
        print(len(table_ukidss), 'rows')

        explore_UKIDSS_LAS(table=table_ukidss)


    if READ_VHS:
        filename_vhs = 'xqz5_xm_VHSv20231101_nn5as.fits'
        filename_vhs = 'xqz5_xm_VHS_VSADR6_nn5as.fits'
        path_vhs = './XMatch/'
        infile = path_vhs + filename_vhs
        print('Reading:', infile)
        table_vhs = Table.read(infile)
        print(len(table_vhs), 'rows')

        explore_VHS(table=table_vhs)


    if READ_VIKING:
        filename_viking = 'xqz5_xm_VIKING_DR5_nn5as.fits'
        path = './XMatch/'
        infile = path + filename_viking
        print('Reading:', infile)
        table_viking = Table.read(infile)
        print(len(table_viking), 'rows')

        explore_VISTA(table=table_viking, viking=True)

        sys.exit()


    lineInfo()
    filename_ls = 'xqz5_xm_ls_dr10_tractor_s_nn5as.fits'
    path_ls = './XMatch/'
    table_ls = read_data(filename=filename_ls, path=path_ls)

    explore_ls(table=table_ls)

    result = xmatch_ls(table_xqz5=table_xqz5,
                       table_ls=table_ls,
                       title_prefix=filename_ls,
                       plotfile_prefix=plotfile_prefix,
                       checkplots=xmatch_checkplots)


    # create html webpage
    mk_webpage(table=result)

    #sys.exit()

    # sanity check plots

    # A4 ratio: 1.414 = (sqrt(2)); 14:10
    # Computer Screen since 2008: 16:9 (1.78:1.0)
    # 4k: 3840 x 2160 1.78
    plt.figure(figsize=(10.0, 6.0))


    # g - r
    # explore differences between no data and -ve fluxes
    # masked means no data; also isnan and is not finite
    # inf

    #
    data = result['mag_g']
    print('mag_g:', len(data))
    print()

    # def explore_masked_data()
    # Only g band has masked values and these are -ve flux values
    # Inf is present when there is no data
    # NaN covers both Inf and Masked

    #help(data)
    print(data.data)
    print(data.mask)

    print(data.description)
    print(data.unit)
    print(data.name)
    print(data.meta)
    print(data.format)
    print()

    imasked_g = ma.is_masked(data)
    itest_is_mask_g = ma.is_mask(data.mask)

    itest_is_masked_g = data.mask
    itest_isnan_g  = np.isnan(data)
    itest_isinf_g  = np.isinf(data)
    itest_isfinite_g  = np.isfinite(data)

    print('Masked values:', ma.count_masked(data), ma.count(data),
          len(data[itest_is_mask_g]), len(data[itest_is_masked_g]))
    print('min, nanmin:',
          np.min(data), np.nanmin(data))
    print('max, nanmax:',
          np.max(data), np.nanmax(data))
    print('Numpy Masked Array min, max:',
          ma.min(data), ma.max(data))
    print('Isfinite min, max:',
          data[np.isfinite(data)].min(),
          data[np.isfinite(data)].max())
    print('isnan:',
          len(data[itest_isnan_g]),
          np.count_nonzero(np.isnan(data)),
          np.count_nonzero(~np.isnan(data)))
    print('isinf:',
          np.count_nonzero(~np.isinf(data)),
          np.count_nonzero(np.isinf(data)))
    print('isfinite:',
          np.count_nonzero(np.isfinite(data)),
          np.count_nonzero(~np.isfinite(data)))

    input('Enter any key to continue... ')

    data = result['mag_r']
    print()

    print(data.data)

    print('Is masked:', ma.is_masked(data))
    if ma.is_masked(data):
        print(data.mask)
    print()

    if ma.is_masked(data):
        imasked_r = ma.is_masked(data)
        itest_is_mask_r = ma.is_mask(data.mask)
        itest_is_masked_r = data.mask

    itest_isnan_r  = np.isnan(data)
    itest_isinf_r  = np.isinf(data)
    itest_isfinite_r  = np.isfinite(data)

    print('mag_r:', len(data))
    print('Masked values:', ma.count_masked(data), ma.count(data))
    print('min, nanmin:',
          np.min(data), np.nanmin(data))
    print('max, nanmax:',
          np.max(data), np.nanmax(data))
    print('Numpy Masked Array min, max:',
          ma.min(data), ma.max(data))
    print('Isfinite min, max:',
          data[np.isfinite(data)].min(),
          data[np.isfinite(data)].max())

    print('isnan:',
          np.count_nonzero(np.isnan(data)),
          np.count_nonzero(~np.isnan(data)))
    print('isinf:',
          np.count_nonzero(~np.isinf(data)),
          np.count_nonzero(np.isinf(data)))
    print('isfinite:',
          np.count_nonzero(np.isfinite(data)),
          np.count_nonzero(~np.isfinite(data)))

    input('Enter any key to continue... ')

    data = result['mag_i']

    print()
    print(data.data)

    print('Is masked:', ma.is_masked(data))
    if ma.is_masked(data):
        print(data.mask)

    print()
    print('mag_i:', len(data))
    print('Masked values:', ma.count_masked(data), ma.count(data))
    print('min, nanmin:',
          np.min(data), np.nanmin(data))
    print('max, nanmax:',
          np.max(data), np.nanmax(data))
    print('Numpy Masked Array min, max:',
          ma.min(data), ma.max(data))
    print('Isfinite min, max:',
          data[np.isfinite(data)].min(),
          data[np.isfinite(data)].max())

    print('isnan:',
          np.count_nonzero(np.isnan(data)),
          np.count_nonzero(~np.isnan(data)))
    print('isinf:',
          np.count_nonzero(~np.isinf(data)),
          np.count_nonzero(np.isinf(data)))
    print('isfinite:',
          np.count_nonzero(np.isfinite(data)),
          np.count_nonzero(~np.isfinite(data)))

    print()
    input('Enter any key to continue... ')

    data = result['mag_z']
    print('mag_z:', len(data))
    print()

    print(data.data)
    print('Is masked:', ma.is_masked(data))
    if ma.is_masked(data):
        print(data.mask)

    print()

    print('Masked values:', ma.count_masked(data), ma.count(data))
    print('min, nanmin:',
          np.min(data), np.nanmin(data))
    print('max, nanmax:',
          np.max(data), np.nanmax(data))
    print('Numpy Masked Array min, max:',
          ma.min(data), ma.max(data))
    print('Isfinite min, max:',
          data[np.isfinite(data)].min(),
          data[np.isfinite(data)].max())

    print('isnan:',
          np.count_nonzero(np.isnan(data)),
          np.count_nonzero(~np.isnan(data)))
    print('isinf:',
          np.count_nonzero(~np.isinf(data)),
          np.count_nonzero(np.isinf(data)))
    print('isfinite:',
          np.count_nonzero(np.isfinite(data)),
          np.count_nonzero(~np.isfinite(data)))

    input('Enter any key to continue... ')


    # g-r analysis

    # 7 objects masked mag g values due to -ve flux
    # 8 objects have no r measurement
    # 67 objects remain

    plt.subplot(121)
    plot_title = 'XQz5_LSSDR10_20240119.fits'
    xdata = result['MgII Redshift']
    ydata = result['mag_g'] - result['mag_r']

    ndata = len(xdata)

    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)

    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.title(plot_title)
    plt.legend(loc='upper left')

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_g - mag_r')

    #plt.show()


    plt.subplot(122)
    xdata = result['mag_r']
    ydata = result['mag_g'] - result['mag_r']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.legend(loc='upper left')
    plt.title(plot_title)
    plt.xlabel('mag_r')
    plt.ylabel('mag_g - mag_r')


    plotfile = 'XQz5_LSSDR10_20240119_mag_r_v_mag_g-mag_r.png'
    lineInfo()
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    plt.show()


    # Redshift versus mag
    plt.figure(figsize=(10.0, 6.0))
    plt.subplot(121)
    plot_title = 'XQz5_LSSDR10_20240119.fits'
    xdata = result['MgII Redshift']
    ydata = result['mag_g']

    ndata = len(xdata)

    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)

    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.title(plot_title)
    plt.legend(loc='upper left')

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_g')

    #plt.show()


    plt.subplot(122)
    xdata = result['MgII Redshift']
    ydata = result['mag_r']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.legend(loc='upper left')
    plt.title(plot_title)

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_g')


    plotfile = 'XQz5_LSSDR10_20240119_Redshift_v_mag_g+mag_r.png'
    lineInfo()
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    plt.show()


    # now look at g-i
    plt.figure(figsize=(10.0, 6.0))

    plt.subplot(121)

    xdata = result['MgII Redshift']
    ydata = result['mag_g'] - result['mag_i']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.title(plot_title)
    plt.legend(loc='upper left')

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_g - mag_i')


    plt.subplot(122)

    xdata = result['mag_i']
    ydata = result['mag_g'] - result['mag_i']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.legend(loc='upper left')
    plt.title(plot_title)
    plt.xlabel('mag_i')
    plt.ylabel('mag_g - mag_i')

    plotfile = 'XQz5_LSSDR10_20240119_mag_i_v_mag_g-mag_i.png'
    lineInfo()
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    plt.show()


    # g - z
    plt.figure(figsize=(10.0, 6.0))
    plt.subplot(121)

    xdata = result['MgII Redshift']
    ydata = result['mag_g'] - result['mag_z']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.title(plot_title)
    plt.legend(loc='upper left')

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_g - mag_z')


    plt.subplot(122)

    xdata = result['mag_z']
    ydata = result['mag_g'] - result['mag_z']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.legend(loc='upper left')
    plt.title(plot_title)
    plt.xlabel('mag_z')
    plt.ylabel('mag_g - mag_z')

    plotfile = 'XQz5_LSSDR10_20240119_mag_z_v_mag_g-mag_z.png'
    lineInfo()
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    plt.show()




    # r - z
    plt.figure(figsize=(10.0, 6.0))
    plt.subplot(121)

    xdata = result['MgII Redshift']
    ydata = result['mag_r'] - result['mag_z']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.title(plot_title)
    plt.legend(loc='upper left')

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_r - mag_z')


    plt.subplot(122)

    xdata = result['mag_z']
    ydata = result['mag_r'] - result['mag_z']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.legend(loc='upper left')
    plt.title(plot_title)
    plt.xlabel('mag_z')
    plt.ylabel('mag_r - mag_z')

    plotfile = 'XQz5_LSSDR10_20240119_mag_z_v_mag_r-mag_z.png'
    lineInfo()
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    plt.show()



    # i - z
    plt.figure(figsize=(10.0, 6.0))
    plt.subplot(121)

    xdata = result['MgII Redshift']
    ydata = result['mag_i'] - result['mag_z']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.title(plot_title)
    plt.legend(loc='upper left')

    plt.xlabel('MgII Redshift')
    plt.ylabel('mag_i - mag_z')


    plt.subplot(122)

    xdata = result['mag_z']
    ydata = result['mag_i'] - result['mag_z']

    ndata = len(xdata)
    ndata_finite = np.count_nonzero(
        (np.isfinite(xdata) & np.isfinite(ydata)))
    label = str(ndata) + '/' + str(ndata_finite)
    print('label:', label)
    plt.plot(xdata, ydata, '.r',
             label=label)

    plt.legend(loc='upper left')
    plt.title(plot_title)
    plt.xlabel('mag_z')
    plt.ylabel('mag_i - mag_z')

    plotfile = 'XQz5_LSSDR10_20240119_mag_z_v_mag_i-mag_z.png'
    lineInfo()
    print('Saving:', plotfile)
    plt.savefig(plotfile)
    logger.info('Saving:' + plotfile)

    plt.show()
