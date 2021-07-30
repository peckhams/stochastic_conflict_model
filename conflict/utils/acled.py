#
#  Copyright (c) 2021, Scott D. Peckham
#  May 2021.  Started on May 17 for Localized Conflict Modeling.
#             First, read GPW population count data as GeoTIFF.
#             Write-up and theoretical results in late May.
#-------------------------------------------------------------------
#
#  read_acled_data()
#
#-------------------------------------------------------------------

import numpy as np
import pandas as pd
import os, os.path

# import time    # (not used yet)

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


