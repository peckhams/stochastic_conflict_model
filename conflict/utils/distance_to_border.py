#
#  Copyright (c) 2022, Scott D. Peckham
#
#  Feb  2022. Wrote a working version (2022-02-08).
#
#-------------------------------------------------------------------

import numpy as np
import scipy.spatial.distance as ssd
import os, os.path, time

from conflict.utils import ncgs_files
from conflict.utils import rtg_files
from conflict.utils import rti_files

#--------------------------------------------------------------------
#
# read_country_code_grid()
# get_distance_to_border()
# save_grid_to_rtg()
#
#--------------------------------------------------------------------
def read_country_code_grid( nc_file=None, REPORT=True):

    #--------------------------------------------------------     
    # Note:  Obtained this grid of country ISO codes from:
    #        https://zenodo.org/record/4457118#.YgLzE33MJUR
    #-------------------------------------------------------- 
    if (nc_file is None):
        nc_dir  = '/Users/peckhams/Dropbox/GitHub/stochastic_conflict_model/data/'
        nc_file = nc_dir + 'countries_gridded_0.1deg_v0.1.nc'
    
    ncgs = ncgs_files.ncgs_file()
    ncgs.open_file( nc_file )
    var_name_list = ncgs.get_var_names()
    if (REPORT):
        print('var_names in netCDF file =' )
        print( var_name_list )

    var_name_list = ncgs.get_var_names( no_dim_vars=True )  ####
    var_index = 0   # (dim vars are now excluded)
    var_name0  = var_name_list[ var_index ]
    print('var_name guess =', var_name0)
    # var_name = 'iso'
    var_name = 'country'
    
    #------------------------------------------    
    # Global grid, with 0.1 degree grid cells
    #------------------------------------------
    # ncols    = 3601
    # nrows    = 1801
    # xres_sec = 360
    # yres_sec = 360
    #--------------------------------------------
    # Use these to set "extent" in plt.imshow()
    #--------------------------------------------
    minlon = -180.05
    maxlon =  180.05
    minlat = -90.05
    maxlat =  90.05
    extent = [minlon, maxlon, minlat, maxlat]
    
    #----------------------------------------------
    # Read grid from nc_file for given time_index
    #----------------------------------------------
    grid = ncgs.ncgs_unit.variables[ var_name ]
    # (nrows, ncols) = grid.shape
    ## print('grid[0,0] =', grid[0,0])
    ## array = np.array( grid )
    array = np.array( grid ).astype('uint8' )  # only 232 unique values
    
    #--------------------
    # Flip the vertical
    #--------------------
    array = np.flipud( array )
     
    if (REPORT):
        print( 'extent = ')
        print( extent )
        print('grid dtype =',  grid.dtype)
        print('grid shape =', grid.shape )
        print('unique values =')
        print( np.unique(grid) )  # 0 through 232

    ncgs.close_file()
    return array
    ## return grid
    
#   read_country_code_grid()
#--------------------------------------------------------------------
def clip_global_grid( grid, minlat, maxlat, minlon, maxlon):

    return grid2

#   clip_global_grid()
#--------------------------------------------------------------------
def get_distance_to_border( stride=10, HORN_OF_AFRICA=True):

    start_time = time.time()
    
    if (HORN_OF_AFRICA):
        CLIP   = True
        minlat = -5.0
        maxlat = 25.0
        minlon = 25.0
        maxlon = 55.0
        s1 = 'Horn_of_Africa'
    else:
        #---------
        # Global
        #---------
        CLIP   = False
        minlat = -90.0
        maxlat =  90.0
        minlon = -180.0
        maxlon = 180.0
        s1 = 'Global'

    codes0 = read_country_code_grid()
    print('original shape =', codes0.shape )
    
    #----------------------------------------------   
    # Use nearest interpolation to make smaller ?
    #----------------------------------------------
    codes = codes0[::stride, ::stride]
    (nrows, ncols) = codes.shape
    if (stride != 1):
        print('reduced shape =', codes.shape )

    #--------------------    
    # Initialize arrays
    #--------------------
    cols = np.arange(ncols)
    rows = nrows - np.arange(nrows)  # (flip y-axis)
    lons = np.linspace(-180.0, 180.0, ncols)   # ncols = 3601
    lats = np.linspace(-90.0, 90.0, nrows)     # nrows = 1801
       
    #-------------------------------------------
    # Clip country code grid to bounding box ?
    #-------------------------------------------
    if (CLIP):
        b1 = np.logical_and(lons >= minlon, lons <= maxlon)
        b2 = np.logical_and(lats >= minlat, lats <= maxlat)
        lons  = lons[ b1 ]
        cols  = cols[ b1 ]
        #--------------------
        lats  = lats[ b2 ]
        rows  = rows[ b2 ]
        #--------------------
        ncols = lons.size
        nrows = lats.size
        print('new ncols, nrows =', ncols, ', ', nrows)
        #-------------------------------------------
        minrow = rows.min()
        maxrow = rows.max()
        mincol = cols.min()
        maxcol = cols.max()
        codes = codes[ minrow:maxrow+1, mincol:maxcol+1 ]  ######
        # codes = codes[ rows, cols ]   # (doesn't work)

    unique_codes = np.unique( codes )
    n_codes = unique_codes.size
    print('Number of country codes found =', n_codes)
    # lon_grid, lat_grid = np.meshgrid( lons, lats )
    xg, yg = np.meshgrid( lons, lats )
    bord_dist = np.zeros( codes.shape, dtype='float32') - 1  #######
    print('codes shape     =', codes.shape )
    print('lon grid shape  =', xg.shape )
    print('lat grid shape  =', yg.shape )
    print('bord_dist shape =', bord_dist.shape )

    #-----------------------------------------    
    # Get new info after reducing & clipping
    #-----------------------------------------
    sec_per_deg = 3600
    xrange   = (maxlon - minlon)
    yrange   = (maxlat - minlat)
    #------------------------------------
    xres_deg = xrange / (ncols-1)  
    xres_sec = xres_deg * sec_per_deg
    #------------------------------------
    yres_deg = yrange / (nrows-1)  
    yres_sec = yres_deg * sec_per_deg
    #------------------------------------
    y_south_edge = minlat - yres_deg/2,
    y_north_edge = maxlat + yres_deg/2,
    x_west_edge  = minlon - xres_deg/2,
    x_east_edge  = maxlon + xres_deg/2,
                                     
    #-------------------------------------------
    # Used: h,edges = np.histogram( bins=233 )
    # to get # of grid cells for each code.
    #-------------------------------------------
    # Skip the code "0" for ocean
    # Skip the code "1" for Antarctica?
    # Skip the code "3" for USSR
    # bord_dist defaults to -1.
    #-------------------------------------------
    big_codes = [0,1,3]   # ocean, Antarctica, USSR
    n_countries = codes.max()
    for code in range(n_countries):   # skip 0 & 1
        if (code not in big_codes):
            print('Working on code =', code)
            w1 = (codes == code)    # boolean array
            w2 = np.invert( w1 )
            pA = np.column_stack( (yg[w1].flatten(), xg[w1].flatten() ) )
            pB = np.column_stack( (yg[w2].flatten(), xg[w2].flatten() ) )
            print('pA.shape =', pA.shape)
            print('pB.shape =', pB.shape)
            dist = ssd.cdist( pA, pB, metric='euclidean' )
            bord_dist[w1] = dist.min( axis=1 )

    #--------------------------------------------------    
    # Save the distance-to-border grid as RTG (FLOAT)
    # Save the country codes grid as RTG (BYTE)
    # IDL "BYTE" type is unsigned (same as "uint").
    #--------------------------------------------------
    save_grid_to_rtg( bord_dist,
                      rtg_file=s1 + '_dist_to_border.rtg',
                      ncols=ncols, nrows=nrows, data_type='FLOAT',
                      xres_sec=xres_sec, yres_sec=yres_sec,
                      y_south_edge=y_south_edge,
                      y_north_edge=y_north_edge,
                      x_west_edge =x_west_edge,
                      x_east_edge =x_east_edge )
    save_grid_to_rtg( codes, rtg_file=s1 + '_country_codes.rtg',
                      ncols=ncols, nrows=nrows, data_type='BYTE',
                      xres_sec=xres_sec, yres_sec=yres_sec,
                      y_south_edge=y_south_edge,
                      y_north_edge=y_north_edge,
                      x_west_edge =x_west_edge,
                      x_east_edge =x_east_edge )

    #----------------------
    # Print final message
    #----------------------
    end_time = time.time()
    run_time = (end_time - start_time)
    run_time = run_time / 60.0
    print('Run time =', run_time, ' [minutes]')

#   get_distance_to_border()
#--------------------------------------------------------------------
def save_grid_to_rtg( grid, rtg_file='TEST.rtg',
                      ncols=3601, nrows=1801, data_type='FLOAT',
                      xres_sec=360, yres_sec=360,
                      y_south_edge=-90.05, y_north_edge=90.05,
                      x_west_edge=-180.05, x_east_edge=180.05):

    print('Current directory =', os.getcwd() )
    print('Saving grid to RTG file...')
    rtg = rtg_files.rtg_file()  
    rti = rti_files.make_info( rtg_file, ncols=ncols, nrows=nrows,
                    xres=xres_sec, yres=yres_sec,
                    data_type=data_type, pixel_geom=0,
                    y_south_edge=y_south_edge,
                    y_north_edge=y_north_edge,
                    x_west_edge=x_west_edge,
                    x_east_edge=x_east_edge,
                    box_units='DEGREES')
    OK = rtg.open_new_file( rtg_file, rti )

    if not(OK):
        print('ERROR during open_new_file().')
        return

    #---------------------------
    # Write a grid to the file
    #---------------------------
    rtg.write_grid( grid, VERBOSE=True )  # calls rtg.close()

#   save_grid_to_rtg()
#--------------------------------------------------------------------




