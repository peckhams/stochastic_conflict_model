#
#  Copyright (c) 2021, Scott D. Peckham
#  May 2021.  Started on May 17 for Localized Conflict Modeling.
#             First, read GPW population count data as GeoTIFF.
#             Write-up and theoretical results in late May.
#-------------------------------------------------------------------
#
#  save_acled_data_to_grid()   # 2021-09-15
#  read_acled_data()           # 2021-05 (original)
#
#-------------------------------------------------------------------

import numpy as np
import pandas as pd
import os, os.path
import time
import rti_files

#-------------------------------------------------------------------
def save_acled_data_to_grid( excel_file='acled.xlsx', data_dir=None,
                     bounds=None, year_range=[1997,2019],
                     dlat_arcsecs=450, dlon_arcsecs=450,
                     # dlat_arcsecs=1800, dlon_arcsecs=1800,
                     HORN_OF_AFRICA=True,
                     rts_file='Horn_of_Africa_deaths_by_month.rts',
                     rtg_file='Horn_of_Africa_deaths_in_month1.rtg'):

    #-----------------------------------------------------
    # Note:  To use this, need:
    #        conda install pandas
    #        conda install openpyxl (for excel files)
    #        DON'T: conda install xlrd  (deprecated)
    #-----------------------------------------------------
    # Note:  ACLED data is for all of Africa.
    #        Year range = 1997 to 2019 (23 years).
    #-----------------------------------------------------
    # Note:  bounds = [ minlon, minlat, maxlon, maxlat ]
    #        For the Horn of Africa:
    #            bounds = [ 25.0, -5.0, 55.0, 25.0]
    #        if (dlat == 1) and (dlon == 1), then
    #            ncols = 30, nrows = 30
    #-----------------------------------------------------
    if (data_dir is None):
        home_dir = os.path.expanduser('~') + os.sep
        data_dir = home_dir + 'Data/ACLED_Data/'
    excel_filepath = data_dir + excel_file

    if (HORN_OF_AFRICA):
        bounds = [25.0, -5.0, 55.0, 25.0]

    #--------------------------------------------
    # Read Excel spreadsheet into pandas object
    # This is the slowest part by far.
    #--------------------------------------------
    start_time = time.time()
    print('Reading ACLED Excel file into pandas...')
    df = pd.read_excel( excel_filepath )
    read_time = (time.time() - start_time)
    print('Time to read Excel file into pandas =', read_time, '[secs].')

    #-------------------------------------------------    
    # Print all column header labels and their types
    #-------------------------------------------------
    # Can also use:  df['YEAR'].dtype
    #-------------------------------------------------  
    # print(df.columns)  # (column header labels)
    # print(df.dtypes)   # (column data type; numpy style)

    #---------------------------------------
    # Extract some columns as numpy arrays
    #-------------------------------------------------------------
    # If actor1 and actor2 are the same, treat as same conflict?
    #-------------------------------------------------------------
    # Use: str(dates[0]) to convert to standard datetime string. 
    #-------------------------------------------------------------       
    event_lats   = df['LATITUDE'].to_numpy()     # (float64)
    event_lons   = df['LONGITUDE'].to_numpy()    # (float64)
    event_deaths = df['FATALITIES'].to_numpy()   # (int64)
    event_dates  = df['EVENT_DATE']              # (datetime64[ns])
    event_years  = df['YEAR'].to_numpy()         # (int64)
    event_actor1 = df['ACTOR1']                  # (str)
    event_actor2 = df['ACTOR2']                  # (str)

    #------------------------------------------------  
    # Create an array of month_nums from EVENT_DATE
    #------------------------------------------------
    datetimes = list(map(str, event_dates))
    ## datetimes = np.char.mod('%0.19s', dates)  # (same; faster?)
    ## datetimes = np.char.mod('%s', dates)  # (full length)
    event_month_nums = np.array([int(x[5:7]) for x in datetimes])
    
    #-----------------------------------------------------
    # Note:  bounds = [ minlon, minlat, maxlon, maxlat ]
    #-----------------------------------------------------
    # Note:  Bounding lats and lons are included in box.
    #-----------------------------------------------------      
    if (bounds is None):
        #---------------------------------------------
        # Find smallest geographic bounding box that
        # contains all of the ACLED events
        #---------------------------------------------
        minlon = np.floor( event_lons.min() )
        maxlon = np.ceil(  event_lons.max() )
        minlat = np.floor( event_lats.min() )
        maxlat = np.ceil(  event_lats.max() )
        IN_BOX = np.ones( event_lats.size, dtype='bool_' )
    else:
        minlon  = bounds[0]
        minlat  = bounds[1]
        maxlon  = bounds[2]
        maxlat  = bounds[3]
        w1 = np.logical_and( event_lats >= minlat, event_lats <= maxlat )
        w2 = np.logical_and( event_lons >= minlon, event_lons <= maxlon )
        IN_BOX = np.logical_and( w1, w2 )

    #------------------------------------------------
    # Convert dlat_arcsecs, dlon_arcsecs to degrees
    #------------------------------------------------
    dlon_deg = (dlon_arcsecs / 3600.)
    dlat_deg = (dlat_arcsecs / 3600.)

    #-----------------------------------------------    
    # Get cols and rows arrays for np.searchsorted
    #-----------------------------------------------
    # Must add 1 for np.linspace, as shown, e.g.
    #    np.linspace(0,2,4+1) = [0,0.5,1,1.5,2]
    #-----------------------------------------------
    nlons = np.int32( (maxlon - minlon) / dlon_deg )
    nlats = np.int32( (maxlat - minlat) / dlat_deg )
    lons  = np.linspace(minlon, maxlon, nlons + 1)
    lats  = np.linspace(minlat, maxlat, nlats + 1)
    ## cols  = np.arange(nlons + 1)  # (not needed)
    ## rows  = np.arange(nlats + 1)

    #----------------------------------------
    # Prepare to write to RTG and RTS files
    #----------------------------------------
    rts_unit = open(rts_file, 'wb')
    rtg_unit = open(rtg_file, 'wb')
 
    #----------------------------   
    # Get array of year numbers
    #----------------------------
    minyear = year_range[0]
    maxyear = year_range[1]
    n_years = (maxyear - minyear) + 1
    year_nums = np.arange(n_years) + minyear
    
    #---------------------------------------------------   
    # Create an RTS file, where each grid in the stack
    # has the number of deaths per cell for one year
    #---------------------------------------------------
    for year in year_nums:
        print('Working on year:', year)
        for month in range(1,13):
            w1 = np.logical_and( IN_BOX, event_years==year )
            w2 = np.logical_and( w1, event_month_nums==month )
            #--------------------------------------------------
            year_month_lons   = event_lons[w2]
            year_month_lats   = event_lats[w2]
            year_month_deaths = event_deaths[w2] 
            n_conflicts       = year_month_lons.size  # (in this year-month)
            
            #------------------------------------------
            # For testing: 73 events have deaths=1369
            #------------------------------------------
#             w3 = (year_month_deaths == 1369)
#             n3 = w3.sum()
#             if (n3 > 0):
#                 print('1369 deaths: year=', year,', month=', month, ', events=', n3)

            #---------------------------------------------------- 
            # Get grid cell row,col for each row in spreadsheet
            # cols must be in [0, nlons-1], inclusive.
            # rows must be in [0, nlats-1], inclusive.
            #--------------------------------------------------------
            # Confirmed with Badme, Eritrea that north-south flip
            # is required, as shown.
            #--------------------------------------------------------
            # Use np.searchsorted, but must use "side='right'" and
            # must substract 1.  For example, if nc = ncols = 4:
            # a = np.linspace(2,4,nc+1) gives [2, 2.5, 3, 3.5, 4]
            #--------------------------------------------------------
            # np.searchsorted(a,2.0) gives 0
            # np.searchsorted(a,2.1) gives 1 (vs. 0)
            # np.searchsorted(a,3.9) gives 4 (vs. 3)
            # np.searchsorted(a,4.0) gives 4 (vs. 3)
            # - If we subtract 1, first line will give -1, but the
            #   min is 0 in this case.  Can fix this with:
            #   np.maximum( cell_rows, 0, cell_rows)  # (in-place)    
            #--------------------------------------------------------
            if (n_conflicts > 0):
                cell_cols = np.searchsorted( lons, year_month_lons) - 1
                cell_rows = np.searchsorted( lats, year_month_lats) - 1
                np.maximum( cell_cols, 0, cell_cols)  # (in-place)
                np.maximum( cell_rows, 0, cell_rows)  # (in-place)
                cell_rows = (nlats - 1 - cell_rows)  # (flip north-south)

            #-------------------------------------------------------------            
            # np.searchsorted(a,2.0, side='right') gives 1 (vs. 0)
            # np.searchsorted(a,2.1, side='right') gives 1 (vs. 0)
            # np.searchsorted(a,3.9, side='right') gives 4 (vs. 3)
            # np.searchsorted(a,4.0, side='right') gives 5 (vs. 3)
            # - If we subtract 1, last line will give 4, but the
            #   max is 3 in this case.  Can fix this with:
            #   np.minimum( cell_rows, nlats-1, cell_rows)  # (in-place)
            #-------------------------------------------------------------
            # if (n_conflicts > 0):
            #     cell_cols = np.searchsorted( lons, year_month_lons, side='right' ) - 1
            #     cell_rows = np.searchsorted( lats, year_month_lats, side='right' ) - 1
            #     np.minimum( cell_cols, nlons-1, cell_cols)  # (in-place)
            #     np.minimum( cell_rows, nlats-1, cell_rows)  # (in-place)
            #     cell_rows = (nlats - 1 - cell_rows)  # (flip north-south)

            #--------------
            # For testing
            #--------------
#             if (n_conflicts > 0):
#                 ## print('cell_cols.size =', cell_cols.size)
#                 ## print('cell_rows.size =', cell_rows.size)
#                 min_cell_col = cell_cols.min()
#                 max_cell_col = cell_cols.max()
#                 min_cell_row = cell_rows.min()
#                 max_cell_row = cell_rows.max()
#                 if (min_cell_col < 0) or (max_cell_col > nlons-1):
#                     print('ERROR:  Some cell cols are out of range.')
#                     print('min_cell_col, max_cell_col =', min_cell_col, max_cell_col)
#                 if (min_cell_row < 0) or (max_cell_row > nlats-1):
#                     print('ERROR:  Some cell rows are out of range.') 
#                     print('min_cell_row, max_cell_row =', min_cell_row, max_cell_row)

            #----------------------------------------------            
            # If there are no deaths in this year-month,
            # we still want to write out a grid of zeros.
            # For loop will be skipped.
            #----------------------------------------------         
            death_grid  = np.zeros( (nlats, nlons), dtype='float32')
            for k in range(n_conflicts):
                row = cell_rows[k]
                col = cell_cols[k]
                n_deaths = year_month_deaths[k]
                ### print('n_deaths =', n_deaths)
                death_grid[ row, col ] += n_deaths
#               # For testing
#                 if (year == 1999) and (month == 2):
#                     print('row, col =', row, col)
#                     print('death_grid[row, col] =', death_grid[row,col])
#                     print()

            #-------------------------------            
            # Write death_grid to RTS file
            #-------------------------------
            death_grid.tofile( rts_unit )

            #----------------------------------        
            # Make an RTG file for first year
            #----------------------------------
            if (year == minyear) and (month == 1):
                death_grid.tofile( rtg_unit)
                rtg_unit.close()
            
    #---------------------       
    # Close the RTS file
    #---------------------
    rts_unit.close()
    print('Finished writing RTS file:')
    print('   ' + rts_file)

    #----------------------------------
    # Create an RTI file for RTG file
    #----------------------------------
    info = rti_files.make_info(
               grid_file=rtg_file, ncols=nlons, nrows=nlats,
               xres=dlon_arcsecs, yres=dlat_arcsecs,
               #---------------------------------------------
               data_source='Stochastic Conflict Model',
               data_type='FLOAT', byte_order='LSB',
               # byte_order=get_rti_byte_order(),
               pixel_geom=0, zres=0.01, z_units='unknown',
               y_south_edge=minlat, y_north_edge=maxlat,
               x_west_edge =minlon, x_east_edge =maxlon,
               box_units='DEGREES')
    rti_files.write_info(rtg_file, info)

    #----------------------------------
    # Create an RTI file for RTS file
    #----------------------------------
    info = rti_files.make_info(
               grid_file=rts_file, ncols=nlons, nrows=nlats,
               xres=dlon_arcsecs, yres=dlat_arcsecs,
               #---------------------------------------------
               data_source='Stochastic Conflict Model',
               data_type='FLOAT', byte_order='LSB',
               # byte_order=get_rti_byte_order(),
               pixel_geom=0, zres=0.01, z_units='unknown',
               y_south_edge=minlat, y_north_edge=maxlat,
               x_west_edge =minlon, x_east_edge =maxlon,
               box_units='DEGREES')
    rti_files.write_info(rts_file, info)
            
    #-----------------------    
    # Print final messages
    #-----------------------
    run_time = (time.time() - start_time)
    print('nlons, nlats =', nlons, nlats)
    print('dlon_arcsecs, dlat_arcsecs =', dlon_arcsecs, dlat_arcsecs)
    print('Time to read Excel file into pandas =', read_time, '[secs].')
    print('Total run time =', run_time, '[secs].')
        
#   save_acled_data_to_grid()
#-------------------------------------------------------------------
#   ORIGINAL CODE:  OBSOLETE SOON.
#-------------------------------------------------------------------
def read_acled_data( excel_file='acled.xlsx', data_dir=None,
                     bounds=None, year_range=[1997,2019],
                     dlat=0.5, dlon=0.5, 
                     rts_file='Horn_of_Africa_deaths_by_year.rts',
                     rtg_file='Horn_of_Africa_deaths_in_year1.rtg'):

    #-----------------------------------------------------
    # Note:  To use this, need:
    #        conda install pandas
    #        conda install openpyxl (for excel files)
    #        DON'T: conda install xlrd  (deprecated)
    #-----------------------------------------------------
    # Note:  Add ability to restrict to bounding box.
    #        ACLED data is for all of Africa, from 1997.
    #-----------------------------------------------------
    # Note:  bounds = [ minlon, minlat, maxlon, maxlat ]
    #        For the Horn of Africa:
    #        bounds = [ 25.0, -5.0, 55.0, 25.0]
    #        if (dlat == 1) and (dlon == 1), then
    #            ncols = 30, nrows = 30
    #-----------------------------------------------------
    if (data_dir is None):
        home_dir = os.path.expanduser('~') + os.sep
        data_dir = home_dir + 'Data/ACLED_Data/'
    excel_filepath = data_dir + excel_file

    df = pd.read_excel( excel_filepath )

    # print(df.columns)  # (print the column header labels)
    # print(df.dtypes)   # (print data type of each column)
    
    lats   = df['LATITUDE'].to_numpy()       # (float64)
    lons   = df['LONGITUDE'].to_numpy()      # (float64)
    deaths = df['FATALITIES'].to_numpy()   # (int64)
    dates  = df['EVENT_DATE']
    years  = df['YEAR'].to_numpy()  # (int64)
  
    #-----------------------------------------------------
    # Note:  bounds = [ minlon, minlat, maxlon, maxlat ]
    #-----------------------------------------------------  
    if (bounds is None):
        minlon = np.floor( lons.min() )
        maxlon = np.ceil( lons.max() )
        minlat = np.floor( lats.min() )
        maxlat = np.ceil( lats.max() )
        IN_BOX = np.ones( lats.size, dtype='bool_' )
    else:
        minlon  = bounds[0]
        minlat  = bounds[1]
        maxlon  = bounds[2]
        maxlat  = bounds[3]
        w1 = np.logical_and( lats >= minlat, lats <= maxlat )
        w2 = np.logical_and( lons >= minlon, lons <= maxlon )
        IN_BOX = np.logical_and( w1, w2 )

    nlons = np.int32( (maxlon - minlon) / dlon )
    nlats = np.int32( (maxlat - minlat) / dlat )
    cols  = np.linspace(minlon, maxlon, nlons)
    rows  = np.linspace(minlat, maxlat, nlats)

    #---------------------------------------------------   
    # Create an RTS file, where each grid in the stack
    # has the number of deaths per cell for one year
    #---------------------------------------------------
    rts_unit = open(rts_file, 'wb')
    rtg_unit = open(rtg_file, 'wb')
    minyear = year_range[0]
    maxyear = year_range[1]
    n_years = (maxyear - minyear) + 1
    year_nums = np.arange(n_years) + minyear
    for year in year_nums:
        print('Working on year:', year)
        w3 = np.logical_and( IN_BOX, years==year )
        year_lons   = lons[w3]
        year_lats   = lats[w3]
        year_deaths = deaths[w3] 
        # year_dates = dates[w3]
        #------------------------------------------------- 
        # Get grid cell row,col for a row in spreadsheet
        #-------------------------------------------------        
        cell_cols   = np.searchsorted( cols, year_lons )
        cell_rows   = np.searchsorted( rows, year_lats )
        death_grid  = np.zeros( (nlats, nlons), dtype='float32')
        n_conflicts = year_lons.size
        for k in range(n_conflicts):
            death_grid[ cell_rows[k], cell_cols[k]] += year_deaths[k]
        # write death_grid to RTS file
        death_grid.tofile( rts_unit )

        #----------------------------------        
        # Make an RTG file for first year
        #----------------------------------
        if (year == minyear):
            death_grid.tofile( rtg_unit)
            rtg_unit.close()
            
    #---------------------       
    # Close the RTS file
    #---------------------
    rts_unit.close()
    print('Finished writing RTS file:')
    print(rts_file)
    return
    #----------------------------------------------------
    # Remaining code is for TOTAL deaths in year_range.
    #----------------------------------------------------
   

    #---------------------------------------------------- 
    # Get grid cell row,col for each row in spreadsheet
    #----------------------------------------------------
    # https://stackoverflow.com/questions/51573164/
    # dividing-geographical-region-into-equal-sized-
    # grid-and-retrieving-indexing-posit
    #----------------------------------------------------
    cell_cols = np.searchsorted( cols, lons )
    cell_rows = np.searchsorted( rows, lats )

    #--------------------------------------------------    
    # Create a geospatial grid of all conflict deaths
    #---------------------------------------------------
    # This for loop is slow, but we need to accumulate
    # the death count from different events.
    # Need a numpy trick to avoid the for loop.
    #-------------------------------------------------------
    # NOTE!  Rows 1 to 74 in the ACLED data for different
    #        months in 1999 all have the same number of
    #        deaths, 1369, which is the max over all rows.
    #        Rows 81 to 198 all have exactly 1000 deaths.
    #        Rows 256-258 are identical except for dates
    #        (3 consecutive days) and each has 400 deaths.
    #        The cause of this overcounting is not clear.
    #-------------------------------------------------------   
    n_conflicts = lons.size
    for k in range(n_conflicts):
        death_grid[ cell_rows[k], cell_cols[k]] += deaths[k]
    return death_grid

    #--------------------------------------------------    
    # Create a geospatial grid of all conflict deaths
    # that occurred in a given year
    #--------------------------------------------------
#     death_grid  = np.zeros( (nlats, nlons), dtype='float64')
#     n_conflicts = lons.size
#     for k in range(n_conflicts):
#         death_grid[ cell_rows[k], cell_cols[k]] += deaths[k]
            
#   read_acled_data()
#-------------------------------------------------------------------


