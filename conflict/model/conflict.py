
#  Copyright (c) 2021, Scott D. Peckham
#
#  Sept 2021. create_rti_file().  Convert RTS to netCDF.
#  Aug  2021. 
#  July 2021. Jupyter notebook with ipywidgets GUI and ability
#             to visualize model output.
#  June 2021. Further development and testing.
#             Prepared as a Python package.
#  May  2021. Started on May 17 for Localized Conflict Modeling.
#             First, read GPW population count data as GeoTIFF.
#             Write-up and theoretical results in late May.
#-------------------------------------------------------------------
#
#  test1()
#  test2()
#
#  class conflict()
#      initialize()
#      read_config_file()
#      create_rti_file()
#      initialize_U()
#      initialize_C1()
#      initialize_C2()
#-----------------------
#      update()
#      update_U()
#      update_C1()
#      update_C2()
#-----------------------
#      update_p()
#      update_S()
#      update_S1()
#      update_S2()       # (obsolete soon?)
#      update_S_old()    # (obsolete soon?)
#      update_time()
#--------------------------------------
#      get_neighbor_cols_and_rows()
#      get_neighbor_values()
#      spread_conflicts_local1()
#      spread_conflicts_nonlocal1()
#      finalize()
#      run_model()
#
#-------------------------------------------------------------------
import numpy as np
import numpy.random as rn

import time
import os, os.path
from conflict.utils import geotiff as geo
from conflict.utils import rti_files

# from conflict.utils import acled as ac

#-------------------------------------------------------------------
def test1():

    cfg_file = 'conflict.cfg'
    c = conflict()
    c.run_model( cfg_file )

#   test1()
#-------------------------------------------------------------------
def test2():

    pop_grid = geo.read_geotiff()

#   test2()
#-------------------------------------------------------------------
def test3( SUBSAMPLE=False ):

    #--------------------------------
    # Use Horn of Africa as a test.
    #--------------------------------
    in_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
    in_file  = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
    in_file  = in_dir + in_file
    # Bounds = [ minlon, minlat, maxlon, maxlat ]
    out_bounds = [ 25.0, -5.0, 55.0, 25.0]
            
    if not(SUBSAMPLE):
        #--------------------------------------
        # 3600 cols x 3600 rows, 30 arcseconds
        #--------------------------------------
        out_file = 'Horn_of_Africa_GPW-v4_pop_count_2020_30sec.tif'
        out_file = in_dir + out_file
        out_xres_sec = None    # will default to in_xres_sec
        out_yres_sec = None    # will default to in_yres_sec
        print('Reading & clipping GeoTIFF file...')
    else:
        #--------------------------------------
        # 360 cols x 360 rows, 300 arcseconds
        #--------------------------------------
        out_file = 'Horn_of_Africa_GPW-v4_pop_count_2020_450sec.tif'
        out_file = in_dir + out_file
        out_xres_sec = 450.0  # (15 times lower resolution)
        out_yres_sec = 450.0  # (15 times lower resolution)
        print('Reading, clipping & subsampling GeoTIFF file...') 
        
        #--------------------------------------
        # 360 cols x 360 rows, 300 arcseconds
        #--------------------------------------
#         out_file = 'Horn_of_Africa_GPW-v4_pop_count_2020_300sec.tif'
#         out_file = in_dir + out_file
#         out_xres_sec = 300.0  # (10 times lower resolution)
#         out_yres_sec = 300.0  # (10 times lower resolution)
#         print('Reading, clipping & subsampling GeoTIFF file...')   

    geo.regrid_geotiff(in_file=in_file, out_file=out_file, 
           out_bounds=out_bounds,
           out_xres_sec=out_xres_sec,
           out_yres_sec=out_yres_sec,
           RESAMPLE_ALGO='bilinear', REPORT=True)

#   test3()                   
#-------------------------------------------------------------------
class conflict():
    #---------------------------------------------------------------
    def initialize( self, cfg_file=None ):

        home_dir = os.path.expanduser('~') + os.sep
                    
        if (cfg_file is not None):
            #-----------------------------------
            # Read params from the config file
            #-----------------------------------
            self.cfg_file = cfg_file
            self.read_config_file()
        else:
            self.n_steps  = 100       
            self.n_cols   = 240
            self.n_rows   = 240       
            #----------------------------------------------
            self.U_file   = ''     # (To use uniform U)            
            self.C1_file  = ''
            self.C2_file  = ''
            #-----------------------------------------------------                    
            # self.in_dir   = home_dir + 'Conflict/Data/GPW-v4/'
            self.in_dir   = home_dir + 'Conflict/Input/'
            self.out_dir  = home_dir + 'Conflict/Output/' 
            cfg_file      = self.in_dir  + 'conflict.cfg'
            self.out_file = self.out_dir + 'conflicts.rts'
            self.IDs_file = self.out_dir + 'conflict_IDs.rts'
            #----------------------------------------------------
            self.c_emerge = 0.001
            ## self.c_emerge = 0.2
            ## self.c_emerge = 0.001  # (must be in (0,1])
            ## self.c_spread = 0.1
            self.c_spread = 0.4
            ## self.c_spread = 0.03
            ## self.c_spread = 0.05
            self.c_spread2 = 0.0
            self.p_resolve = 0.4
            #----------------------------------------
            self.spread_method = 1
            self.time_lag = 1    # (not used yet)
            self.REPORT   = True
    
        self.home_dir     = home_dir
        self.time_index   = 0
        self.n_conflict_cells = 0
        self.grid_shape  = (self.n_rows, self.n_cols)
        ## self.start_time  = time.time()
        ### self.start_ID    = 1
        self.start_index = 0

        #----------------------------------
        # Georeferencing:  Horn of Africa
        #--------------------------------------
        # Move some of this into the CFG file
        #--------------------------------------
        self.min_lat = -5.0
        self.max_lat = 25.0
        self.min_lon = 25.0
        self.max_lon = 55.0
        lon_range_deg = (self.max_lon - self.min_lon)
        lat_range_deg = (self.max_lat - self.min_lat)
        #-----------------------------------------------
        # If n_cols = n_rows = 240, xres = yres = 450.
        # If n_cols = n_rows = 360, xres = yres = 300.
        #-----------------------------------------------        
        self.xres_arcsecs = lon_range_deg * 3600 / self.n_cols
        self.yres_arcsecs = lat_range_deg * 3600 / self.n_rows
        self.create_rti_file()

        #----------------------------
        # Change to input directory
        #----------------------------
        ## os.chdir( self.in_dir )
     
        #-----------------------------
        # Open output files to write
        #-----------------------------
        self.out_unit = open( self.out_file, 'wb')
        self.IDs_unit = open( self.IDs_file, 'wb')

        #--------------------------------------       
        # Make grids with col and row numbers
        #--------------------------------------
        cols   = np.arange( self.n_cols )
        rows   = np.arange( self.n_rows )
        cg, rg = np.meshgrid( cols, rows )
        self.col_grid = cg
        self.row_grid = rg
            
        self.initialize_U()
        self.initialize_C1()
        self.initialize_C2()
        
        #------------------------------------------------------    
        # Initialize to no conflicts
        # S will later contain 1s in grid cells with conflict
        # Initialize durations to zero also.
        # IDs will contain a unique ID for each conflict.
        # Using 'float32' for IDs now for viewing the RTS.
        #------------------------------------------------------
        self.S    = np.zeros( self.grid_shape, dtype='uint8' )
        self.durs = np.zeros( self.grid_shape, dtype='uint32')
        self.IDs  = np.zeros( self.grid_shape, dtype='float32')
        
        #----------------------------------------------------------
        # Create a set of random integer IDs, without replacement
        # so when we colorize, it will look better.
        # (I think this works as an iterator.)
        #----------------------------------------------------------
        max_ID = 10000000
        self.ran_IDs = rn.choice( max_ID, max_ID, replace=False)
        # This next method used built-in random & and problems.
        ### self.ran_IDs = rn.sample( range(1000000), 500000)
        
        #------------------------------------------------
        # Assign ocean grid cells their own, nonzero ID
        # Use new "LAND_SEA_BACKDROP" option (via w1).
        #------------------------------------------------
        w1 = (self.U == 0)   # boolean array
        self.IDs[ w1 ] = max_ID
        # self.IDs[ w1 ] = self.ran_IDs[0]

        self.start_time  = time.time()
        
    #   initialize()
    #---------------------------------------------------------------
    def read_config_file( self, delim="=" ):

        VERBOSE = False   # Set to True for testing
        cfg_file2 = os.path.expanduser( self.cfg_file )
        cfg_unit = open( cfg_file2, 'r')
        while (True):
            line = cfg_unit.readline()
            if (line == ''): break

            #--------------------------------------
            # Extract the variable name or label,
            # which may contain blank spaces
            #--------------------------------------------------
            # strip() removes leading and trailing whitespace
            #--------------------------------------------------
            p = line.find( delim )
            if (p != -1) and not(line.startswith('#')): 
                key   = line[:p].strip()
                value = line[p + 1:]
                value = value.strip()
                #-----------------------------------
                # Type can be string, float or int
                #-----------------------------------
                if (value[0] in ["'", '"']):
                    value = str(value[1:-1])  # strip quotes
                    value = os.path.expanduser( value )
                else:
                    if ('.' in value):
                        value = np.float32(value)
                    else:
                        value = np.int32(value)
                setattr( self, key, value )
                if (VERBOSE):
                    print('key =', key, ', value =', value)
                    print('type(value) =', type(value)) 
            else:
                # Skip over comments & invalid lines
                pass

        cfg_unit.close()
        if (VERBOSE):
            print('Finished reading cfg file:')
            print('  ' + self.cfg_file)
            print()
        
    #   read_config_file()
    #---------------------------------------------------------------    
    def create_rti_file( self ):
    
        rts_file = self.IDs_file
        info = rti_files.make_info(
                  grid_file=rts_file,
                  ncols=self.n_cols, nrows=self.n_rows,
                  xres=self.xres_arcsecs,
                  yres=self.yres_arcsecs,
                  #---------------------------------
                  data_source='Stochastic Conflict Model',
                  data_type='FLOAT',
                  # byte_order='LSB',
                  # byte_order=get_rti_byte_order(),
                  pixel_geom=0,
                  zres=0.01, z_units='unknown',
                  y_south_edge=self.min_lat,
                  y_north_edge=self.max_lat,
                  x_west_edge =self.min_lon,
                  x_east_edge =self.max_lon,
                  box_units='DEGREES')
 
        rti_files.write_info(rts_file, info)
 
    #   create_rti_file()
    #---------------------------------------------------------------
    def initialize_U( self ):

        #-----------------------------------    
        # Start with U = a population grid
        #-----------------------------------
        if (self.U_file != ''):
            self.U = geo.read_geotiff(in_file=self.U_file,
                                      REPORT=True)
            # In case of negative nodata value
            np.maximum(self.U, 0.0, self.U)   # (in place)
        else:        
            #---------------------    
            # Use a grid of ones
            #---------------------
            self.U = np.ones( self.grid_shape, dtype='float32' )
 
        #-----------------------------------
        # Disallow conflict on the 4 edges
        #-----------------------------------
        self.U[0,:]               = 0.0
        self.U[self.n_rows - 1,:] = 0.0
        self.U[:,0]               = 0.0
        self.U[:,self.n_cols - 1] = 0.0
                      
    #   initialize_U()
    #---------------------------------------------------------------
    def initialize_C1( self ):

        if (self.C1_file != ''):
            self.C1 = geo.read_geotiff(in_file=self.C1_file,
                                       REPORT=True)
            # In case of negative nodata value
            np.maximum(self.C1, 0.0, self.C1)   # (in place)
        else:        
            #---------------------    
            # Use a grid of ones
            #---------------------
            self.C1 = np.ones( self.grid_shape, dtype='float32' )

        #---------------------------------------------
        # Disallow conflict spreading on the 4 edges
        #---------------------------------------------
        self.C1[0,:]               = 0.0
        self.C1[self.n_rows - 1,:] = 0.0
        self.C1[:,0]               = 0.0
        self.C1[:,self.n_cols - 1] = 0.0
                
    #   initialize_C1()
    #---------------------------------------------------------------
    def initialize_C2( self ):
    
        if (self.C2_file != ''):
            self.C2 = geo.read_geotiff(in_file=self.C2_file,
                                       REPORT=True)
            # In case of negative nodata value
            np.maximum(self.C2, 0.0, self.C2)   # (in place)
        else:        
            #---------------------    
            # Use a grid of ones
            #---------------------
            self.C2 = np.ones( self.grid_shape, dtype='float32' )

        #---------------------------------------------
        # Disallow conflict spreading on the 4 edges
        #---------------------------------------------
        self.C2[0,:]               = 0.0
        self.C2[self.n_rows - 1,:] = 0.0
        self.C2[:,0]               = 0.0
        self.C2[:,self.n_cols - 1] = 0.0
        
    #   initialize_C2()
    #---------------------------------------------------------------
    def update( self ):

        self.update_U()
        self.update_C1()
        self.update_C2()
        #-------------------  
        self.update_p()
        self.update_S()      # also updates IDs
        ## self.update_S1()  # same speed as update_S()
        ## self.update_S2()
        self.update_time()

    #   update()
    #---------------------------------------------------------------
    def update_U( self ):
    
        pass

    #   update_U()
    #---------------------------------------------------------------
    def update_C1( self ):
    
        pass

    #   update_C1()
    #---------------------------------------------------------------
    def update_C2( self ):
    
        pass

    #   update_C2()    
    #---------------------------------------------------------------
    def update_p( self ):

        #----------------------------------------------------------
        # Note:  p is the probability that a conflict emerges in
        #        a grid cell, and is a function of the unrest, U.
        #        In order for p to be a probability, in (0,1],
        #        we need 0 < c_emerge <= 1.
        #----------------------------------------------------------
        self.p_emerge = (self.c_emerge / self.U.max()) * self.U
        
    #   update_p()
    #---------------------------------------------------------------
    def update_S( self ):

        #-----------------------------------------------------------
        # Note:  The previous version of this method generated
        #        Geometric random variables to model conflict
        #        durations.  This new version does not track
        #        durations explicitly, but should produce exactly
        #        the same result.  Here, any conflict ends in the
        #        kth time interval with fixed probability, p.
        #        This is modeled with a Bernoulli random variable.
        #-----------------------------------------------------------
        # Note:  A Bernoulli random variable takes the value 1
        #        with probability p and 0 with probability (1-p).
        #        It is a special case of a Binomial r.v. (n=1).
        #        np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 
                
        #--------------------------------------------------
        # Initiate new conflicts in cells with no conflict
        #-----------------------------------------------------
        # Generate Bernoulli random variables with parameter
        # p_emerge, and initiate conflicts where B1=1. 
        #-----------------------------------------------------  
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #-----------------------------------------------------
        B1 = np.random.binomial(1, self.p_emerge)
        w2 = np.logical_and( self.S == 0, B1 == 1 )
        n2 = w2.sum()

        #------------------------------------------------
        # Resolve some conflicts in cells with conflict
        #-----------------------------------------------------        
        # Generate Bernoulli random variables with parameter
        # p_resolve and terminate conflicts that get B2=1.
        # Conflict durations will then turn out to be
        # Geometric random variables, same parameter.
        #-----------------------------------------------------
        B2 = np.random.binomial(1, self.p_resolve, size=self.grid_shape)
        w3 = np.logical_and( self.S == 1, B2 == 1 )
        n3 = w3.sum()        
 
        #------------------------------------       
        # Perform the required updates to S
        #------------------------------------
        i = self.start_index
        self.S[ w2 ]   = B1[ w2 ]
        self.IDs[ w2 ] = self.ran_IDs[i:i + n2]
        self.start_index += n2
        self.n_conflict_cells += n2
        #---------------------------------------------
        self.S[ w3 ] = (1 - B2[ w3 ])
        self.n_conflict_cells -= n3        
         # Reset IDs to zero where resolved (i.e. B2 = 1).
        self.IDs[ w3 ] = self.IDs[w3] * (1 - B2[w3])        
        #---------------------------------------------       
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of new conflict cells =', n2)
            print('Number of resolved conflicts =', n3)
            if (self.spread_method == 0):   
                print()
 
        #------------------------------------------   
        # Attempt to spread the conflicts locally
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts_local1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass

        #---------------------------------------------   
        # Attempt to spread the conflicts nonlocally
        #---------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts_nonlocal1()
                 
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S()
    #---------------------------------------------------------------
    def update_S1( self ):

        #-----------------------------------------------------------
        # Note:  The previous version of this method generated
        #        Geometric random variables to model conflict
        #        durations.  This new version does not track
        #        durations explicitly, but should produce exactly
        #        the same result.  Here, any conflict ends in the
        #        kth time interval with fixed probability, p.
        #        This is modeled with a Bernoulli random variable.
        #-----------------------------------------------------------
        # Note:  A Bernoulli random variable takes the value 1
        #        with probability p and 0 with probability (1-p).
        #        It is a special case of a Binomial r.v. (n=1).
        #        np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 

        #----------------------------------------------
        # Make a copy of self.S, i.e. S(k) vs. S(k+1)
        #----------------------------------------------
        # Otherwise, conflicts may be resolved in the
        # same time step as they were initiated.
        # Could also apply self.S updates at end.
        #----------------------------------------------
        S = self.S.copy()   # This may be costly
                
        #--------------------------------------------------
        # Initiate new conflicts in cells with no conflict
        #-----------------------------------------------------
        # Generate Bernoulli random variables with parameter
        # p_emerge, and initiate conflicts where B1=1. 
        #-----------------------------------------------------  
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #-----------------------------------------------------
        B1 = np.random.binomial(1, self.p_emerge)
        w2 = np.logical_and( S == 0, B1 == 1 )
        n2 = w2.sum()
        i = self.start_index
        self.S[ w2 ]   = B1[ w2 ]
        self.IDs[ w2 ] = self.ran_IDs[i:i + n2]
        self.start_index += n2
        self.n_conflict_cells += n2

        #------------------------------------------------
        # Resolve some conflicts in cells with conflict
        #-----------------------------------------------------        
        # Generate Bernoulli random variables with parameter
        # p_resolve and terminate conflicts that get B2=1.
        # Conflict durations will then turn out to be
        # Geometric random variables, same parameter.
        #-----------------------------------------------------
        B2 = np.random.binomial(1, self.p_resolve, size=self.grid_shape)
        w3 = np.logical_and( S == 1, B2 == 1 )
        n3 = w3.sum()
        self.S[ w3 ] = (1 - B2[ w3 ])
        self.n_conflict_cells -= n3        
         # Reset IDs to zero where resolved (i.e. B2 = 1).
        self.IDs[ w3 ] = self.IDs[w3] * (1 - B2[w3])
 
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of new conflicts =', n2)             
            print('Number of resolved conflicts =', n3)   
            print()

        #------------------------------------------   
        # Attempt to spread the conflicts locally
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts_local1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass

        #---------------------------------------------   
        # Attempt to spread the conflicts nonlocally
        #---------------------------------------------
        # if (self.spread_method == 1):
        #     self.spread_conflicts_nonlocal1()
                             
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S1()
    #---------------------------------------------------------------
    def update_S2( self ):

        #-----------------------------------------------------------
        # Note:  The previous version of this method generated
        #        Geometric random variables to model conflict
        #        durations.  This new version does not track
        #        durations explicitly, but should produce exactly
        #        the same result.  Here, any conflict ends in the
        #        kth time interval with fixed probability, p.
        #        This is modeled with a Bernoulli random variable.
        #-----------------------------------------------------------
        # Note:  A Bernoulli random variable takes the value 1
        #        with probability p and 0 with probability (1-p).
        #        It is a special case of a Binomial r.v. (n=1).
        #        np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 

        #---------------------------------------         
        # Find cells with and without conflict
        #---------------------------------------
        w1 = (self.S == 1)
        w0 = np.invert( w1 )
        n1 = w1.sum()
        n0 = w0.sum()
           
        ## S = self.S.copy()  ########
             
        #--------------------------------------------------
        # Initiate new conflicts in cells with no conflict
        #--------------------------------------------------    
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #--------------------------------------------------
        B1 = np.random.binomial(1, self.p_emerge[w0])
        w2 = (B1 == 1)
        n2 = w2.sum()
        i = self.start_index
        self.S.flat[ w2 ]   = B1[w2]
        self.IDs.flat[ w2 ] = self.ran_IDs[i:i + n2]
        self.start_index += n2
        self.n_conflict_cells += n2
       
        if (self.REPORT):
            print('Number of new conflicts =', n2)

        #------------------------------   
        # Update S with new conflicts
        #------------------------------
#         if (n3 > 0):
#             i = self.start_index
#             self.S[ w3 ]   = 1          # (for method 1)
#             self.IDs[ w3 ] = self.ran_IDs[i:i + n3]
#             #---------------------------------------------
#             #### self.S[ w0 ]   = B1   # (for method 2)
#             #### self.IDs[ w0 ] = self.ran_IDs[i:i + n3]
#             #---------------------------------------------
#             self.start_index += n3
#             self.n_conflict_cells += n3

        #------------------------------------------------
        # Resolve some conflicts in cells with conflict
        #-----------------------------------------------------        
        # Generate Bernoulli random variables, and terminate
        # conflicts that get B2=1.  Conflict durations will
        # then turn out to be Geometric random variables.
        #-----------------------------------------------------
        B2 = np.random.binomial(1, self.p_resolve, size=n1)
        w3 = (B2 == 1)
        n3 = w3.sum()
        self.S.flat[ w3 ] = (1 - B2[ w3 ])
        self.n_conflict_cells -= n3        
         # Reset IDs to zero where resolved (i.e. B2 = 1).
        self.IDs.flat[ w3 ] *= (1 - B2[w3])
 
        #-----------------------------------    
        # Update S with resolved conflicts
        #-----------------------------------                
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of resolved conflicts =', n3)   
            
        #----------------------------------   
        # Attempt to spread the conflicts
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts_local1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass
         
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S2()
    #---------------------------------------------------------------
    def update_S_old( self ):

        #-----------------------------------------------------------
        # Notes:  A Bernoulli random variable takes the value 1
        #         with probability p and 0 with probability (1-p).
        #         It is a special case of a Binomial r.v. (n=1).
        #         np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 

        #------------------------------------------         
        # Reduce the existing durations by 1
        # Durations are integers (# of timesteps)
        #------------------------------------------
        w1 = (self.S == 1)
        # n1 = w1.sum()
        self.durs[ w1 ] -= 1

        #---------------------------------------------
        # Have any conflicts reached their duration?
        # If so, set their S value to 0.
        #---------------------------------------------
        # S[w1][w2] works for retrieving values, but
        # it doesn't work for assignments.
        #---------------------------------------------
#         w2 = (self.durs[ w1 ] == 0)  # can be empty
#         self.S[ w1 ][ w2 ] = 0  # conflict is over


        #-----------------------------------------------------
        # METHOD 1:  Works, but seems no faster than METHOD 2
        #-----------------------------------------------------
#         w2 = (self.durs[ w1 ] == 0)  # can be empty
#         self.S[ w1 ] = (1 - w2.astype('uint8'))
        #----------------------------------------------------  
        # METHOD 2
        #----------------------------------------------------        
        w2 = np.logical_and( (self.S == 1),(self.durs == 0) )
        self.S[ w2 ]   = 0
        self.IDs[ w2 ] = 0   # reset the IDs
        n2 = w2.sum()
        self.n_conflict_cells -= n2
                        
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of resolved conflicts =', n2)

        #--------------------------------------------------
        # Initiate new conflicts;  inherit size from p
        #--------------------------------------------------    
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #--------------------------------------------------
        dS = np.random.binomial(1, self.p_emerge)
        dS = dS.astype('uint8')
        w3 = np.logical_and(self.S == 0, dS == 1)
        #-------------------------------------------- 
        # This would allow existing conflicts to be
        # "reseeded" and causes error in counting.
        #--------------------------------------------    
        ### w3 = (dS == 1)
        n3 = w3.sum()
        if (self.REPORT):
            print('Number of new conflicts =', n3)

        if (n3 > 0):
            #------------------------------   
            # Update S with new conflicts
            #------------------------------
            self.S[ w3 ]   = 1
            ## self.IDs[ w3 ] = np.arange( n3 ) + self.start_ID
            ## self.start_ID += n3
            i = self.start_index
            self.IDs[ w3 ] = self.ran_IDs[i:i + n3]
            self.start_index += n3

            ### np.maximum( self.S, dS, self.S)    # in place
            #------------------------------------------         
            # New durations are Geometric random vars
            #------------------------------------------           
            g  = np.random.geometric( self.p_geom, size=n3 )
            self.durs[ w3 ] = g
            self.n_conflict_cells += n3

        #----------------------------------   
        # Attempt to spread the conflicts
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts_local1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass
         
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S_old()
    #---------------------------------------------------------------
    def update_time( self ):

        self.time_index += 1
        
    #   update_time()
    #---------------------------------------------------------------    
    def get_neighbor_cols_and_rows( self, w1, n1, NO_EDGES=False ):

        cols = self.col_grid[ w1 ]
        rows = self.row_grid[ w1 ]

        #------------------------------------------        
        # Exclude edges using in-place operations
        #------------------------------------------
        if (NO_EDGES):
           np.minimum( np.maximum(cols,1,cols), self.n_cols-2, cols)
           np.minimum( np.maximum(rows,1,rows), self.n_rows-2, rows)
        
        #--------------------------------------------------
        # 1st index is over grid cells that have conflict.
        # 2nd index is over the 8 nearest neighbors.
        #--------------------------------------------------
        cn = np.zeros( (n1, 8), dtype='int32')
        cn[:,0] = cols-1
        cn[:,1] = cols
        cn[:,2] = cols+1
        cn[:,3] = cols-1
        cn[:,4] = cols+1
        cn[:,5] = cols-1
        cn[:,6] = cols
        cn[:,7] = cols+1
        #--------------------------------------- 
        rn = np.zeros( (n1, 8), dtype='int32')
        rn[:,0] = rows-1
        rn[:,1] = rows-1
        rn[:,2] = rows-1
        rn[:,3] = rows
        rn[:,4] = rows
        rn[:,5] = rows+1
        rn[:,6] = rows+1 
        rn[:,7] = rows+1 
        #------------------
        self.cn = cn
        self.rn = rn
                    
    #   get_neighbor_cols_and_rows()
    #---------------------------------------------------------------
    def get_neighbor_values( self, var, n1 ):
    
        #----------------------------------------
        # Get values of 8 nearest neighbors
        # vals[k,:] = neighbor values of cell k
        #----------------------------------------
        cn = self.cn
        rn = self.rn
        
        vals = np.zeros( (n1, 8), dtype='float32')
        vals[:,0] = var[rn[:,0], cn[:,0]]    # (top left)
        vals[:,1] = var[rn[:,1], cn[:,1]]    # (top center)
        vals[:,2] = var[rn[:,2], cn[:,2]]    # (top right)
        vals[:,3] = var[rn[:,3], cn[:,3]]    # (left center)
        vals[:,4] = var[rn[:,4], cn[:,4]]    # (right center)
        vals[:,5] = var[rn[:,5], cn[:,5]]    # (bottom left)
        vals[:,6] = var[rn[:,6], cn[:,6]]    # (bottom center)
        vals[:,7] = var[rn[:,7], cn[:,7]]    # (bottom right)
        # vals[:,8] = var[rn[:,8], cn[:,8]]  # (center)
        
        return vals
        
    #   get_neighbor_values()    
    #---------------------------------------------------------------
    def spread_conflicts_local1( self, USE_LOOP=False ):
 
        if (self.c_spread == 0):
            if (self.REPORT):
                print('No local spreading: c_spread = 0.')
                print()
            return

        #-------------------------------------------------   
        # Note:  Can only spread to cells that have S=0.
        #-------------------------------------------------
        w1 = (self.S == 1)
        n1 = w1.sum()
        if (n1 == 0):
            print('No conflicts to spread at time:', self.time_index)
            return

        if (USE_LOOP):
            ID_vals = self.IDs[ w1 ]  #(for the for loop version)
        else:
            ID_vals = np.tile( np.array([self.IDs[w1]]).transpose(), (1,8))
            #------------------
            # This also works
            #------------------
    #         ID_vals = np.zeros((n1,8), dtype='int64')
    #         ID_vals[:,0] = self.IDs[w1]
    #         ID_vals[:,1] = self.IDs[w1]
    #         ID_vals[:,2] = self.IDs[w1]
    #         ID_vals[:,3] = self.IDs[w1]
    #         ID_vals[:,4] = self.IDs[w1]
    #         ID_vals[:,5] = self.IDs[w1] 
    #         ID_vals[:,6] = self.IDs[w1]
    #         ID_vals[:,7] = self.IDs[w1]
      
        #---------------------------------------------
        # Get nearest neighbor values for U, S, & C1
        #---------------------------------------------
        self.get_neighbor_cols_and_rows( w1, n1, NO_EDGES=True)

        #---------------------------------------        
        # One of more of these is needed below
        #---------------------------------------
        ## Sn  = self.get_neighbor_values( self.S,  n1) 
        ## Un  = self.get_neighbor_values( self.U,  n1)
        C1n = self.get_neighbor_values( self.C1, n1)

        #------------------------------------------------        
        # Compute probability of spreading to neighbors
        #------------------------------------------------
        # The "None trick" shown here allows us to do
        # the following for all k at once:
        # pn[k,:]  = C1n[k,:] * (c2 / C1n[k,:].max() )
        # Need c2 = c_spread to be in (0,1].
        # np.amax lets us take the max along an axis.
        #------------------------------------------------
        # NOTE: Un, C1n and pn have shape = (n1, 8)
        # NOTE: pn is initialized & defaults to 0.
        #------------------------------------------------    
        C1n_max = np.amax( C1n, axis=1 )  # a 1D array
        wg = (C1n_max > 0)
        pn = np.zeros(C1n.shape, dtype='float32')
        pn[ wg,: ] = self.c_spread * C1n[wg,:] / (C1n_max[wg,None])

        #------------------------------------
        # Alternate method that uses only U
        #------------------------------------
        # Un_max = np.amax( Un, axis=1 )  # a 1D array
        # wg = (Un_max > 0)
        # pn = np.zeros(Un.shape, dtype='float32')
        # pn[ wg,: ] = self.c_spread * Un[wg,:] / (Un_max[wg,None])
               
        #--------------------------------------
        # Alternate method that uses U and C1
        #--------------------------------------
        # Rn = Un * C1n 
        # Rn_max = np.amax( Rn, axis=1 )  # a 1D array
        # wg = (Rn_max > 0)
        # pn = np.zeros(Rn.shape, dtype='float32')
        # pn[ wg,: ] = self.c_spread * Rn[wg,:] / (Rn_max[wg,None])

        #---------------------------------------------
        # Use Bernoulli r.v.s to determine spreading
        #---------------------------------------------
        cn = self.cn
        rn = self.rn
        n_start = self.n_conflict_cells
        
        if (USE_LOOP):        
            for k in range(n1):
                B  = np.random.binomial(1, pn[k,:])  # (8 r.v.s)
                ## B  = B.astype( 'uint8' )  ######
                w2 = np.logical_and( (self.S[rn[k,:], cn[k,:]] == 0),
                                     (B == 1) )
                n2 = w2.sum()
                #-----------------------------------------
                # Spread conflict to some neighbor cells
                #-----------------------------------------
                # w2 is boolean array, but this is okay
                #----------------------------------------------------
                # Set duration to be a geometric r.v.
                ## g1 = np.random.geometric( self.p_geom, size=n2 )
                ## self.durs[ rn[k,w2], cn[k,w2] ] = g1
                #----------------------------------------------------
                self.S[    rn[k,w2], cn[k,w2] ] = 1
                self.IDs[ rn[k,w2], cn[k,w2] ]  = ID_vals[k]
                self.n_conflict_cells += n2
        else:
            #-------------------------------------
            # Spread conflict without a for loop
            # Much faster, and results look similar
            # See notes below re: overcounting.
            #-------------------------------------
            B  = np.random.binomial(1, pn)  # (8 r.v.s)
            ## B  = B.astype( 'uint8' )  ######
            w2 = np.logical_and( (self.S[rn, cn] == 0), (B == 1) )
            n2 = w2.sum()
            #-----------------------------------------
            # Spread conflict to some neighbor cells
            #-----------------------------------------
            # w2 is boolean array, but this is okay
            #-----------------------------------------
            self.S[ rn[w2], cn[w2] ]   = 1
            self.IDs[ rn[w2], cn[w2] ] = ID_vals[w2]
            #--------------------------------------------
            # Without the for loop, several cells can 
            # spread conflict to the same cell and this
            # next line results in over-counting:
            #    self.n_conflict_cells += n2
            #--------------------------------------------
            w_new = (self.S == 1)
            n_new = (w_new.sum() - n1)
            self.n_conflict_cells += n_new
                     
        if (self.REPORT):
            n_spread = (self.n_conflict_cells - n_start)
            print('Number of spread conflicts =', n_spread)
            print()
        
    #   spread_conflicts_local1()
    #---------------------------------------------------------------
    def spread_conflicts_nonlocal1( self, USE_LOOP=True ):
 
        if (self.c_spread2 == 0):
            if (self.REPORT):
                print('No nonlocal spreading: c_spread2 = 0.')
                print()
            return

        #-------------------------------------------------   
        # Note:  Can only spread to cells that have S=0.
        #-------------------------------------------------
        w1 = (self.S == 1)
        n1 = w1.sum()
        if (n1 == 0):
            print('No conflicts to spread at time:', self.time_index)
            return

        if (USE_LOOP):
            # Usually, only a few of these are unique
            ## ID_vals = self.IDs[ w1 ]  #(for the for loop version)
            
            ID_vals, ID_counts = np.unique( self.IDs[ w1 ], return_counts=True)
            n_IDs = ID_vals.size
#         else:
#             ID_vals = np.tile( np.array([self.IDs[w1]]).transpose(), (1,8))

        #--------------------------------------------------------        
        # Compute probability that a given conflict (many cells
        # with same ID) will spread to any other grid cell
        #--------------------------------------------------------
        # The "size" of a given conflict can be measured as
        # the number of conflict grid cells with the same ID.
        # We could take the probability that a given conflict
        # spreads to remote cells as proportional to this size,
        # because there is more "coverage" of larger conflicts.
        #-------------------------------------------------------
        n_start = self.n_conflict_cells
        ID_count_max = np.max( ID_counts )
        ## product_max = np.max( self.U ) * ID_count_max
                
        for k in range(n_IDs):
            # Use U to exclude ocean & lake grid cells, etc.
            w2 = np.logical_and(self.S == 0, self.U != 0)
            ## w2 = (self.S == 0)   # (no-conflict cell locations)
            n2 = w2.sum()

            #----------------------------------------------            
            # Generate grid of Bernoulli random variables
            # using a grid of probabilities
            #----------------------------------------------
            if (n2 > 0):
                p = np.zeros(self.C2.shape, dtype='float32')
                #---------------------------------------------------
                # Option 1.  Both U and size of conflict affect p.
                #---------------------------------------------------
                product_max = np.max( self.U[w2] ) * ID_count_max
                R = (self.U[w2] * ID_counts[k]) / product_max
                p[ w2 ] = self.c_spread2 * R
                #---------------------------------------------------
                # R = (self.U * ID_counts[k]) / product_max
                # p[ w2 ] = self.c_spread2 * R[ w2 ]
                                                    
                #----------------------------------------------
                # Option 2.  Only size of conflict affects p.
                #----------------------------------------------
                # p[ w2 ] = self.c_spread2 * ID_counts[k] / ID_count_max
                
                B  = np.random.binomial(1, p)  # entire grid
                ## B  = B.astype( 'uint8' )  ######
                w3 = (B == 1)   # (boolean array)
                n3 = w3.sum()

                self.S[ w3 ] = 1
                self.IDs[ w3 ]  = ID_vals[k]
                self.n_conflict_cells += n3
                     
        if (self.REPORT):
            n_spread = (self.n_conflict_cells - n_start)
            print('Number of spread conflicts =', n_spread)
            print()
        
    #   spread_conflicts_nonlocal1()
    #---------------------------------------------------------------
    def spread_conflicts_nonlocal2( self, USE_LOOP=True ):
 
        if (self.c_spread2 == 0):
            if (self.REPORT):
                print('Number of spread conflicts = 0.')
                print('  c_spread2 = 0.')
                print()
            return

        #-------------------------------------------------   
        # Note:  Can only spread to cells that have S=0.
        #-------------------------------------------------
        w1 = (self.S == 1)
        n1 = w1.sum()
        if (n1 == 0):
            print('No conflicts to spread at time:', self.time_index)
            return

        if (USE_LOOP):
            # Usually, only a few of these are unique
            ## ID_vals = self.IDs[ w1 ]  #(for the for loop version)
            
            ID_vals, ID_counts = np.unique( self.IDs[ w1 ], return_counts=True)
            n_IDs = ID_vals.size
#         else:
#             ID_vals = np.tile( np.array([self.IDs[w1]]).transpose(), (1,8))

        #--------------------------------------------------------        
        # Compute probability that a given conflict (many cells
        # with same ID) will spread to any other grid cell
        #--------------------------------------------------------
        # The "size" of a given conflict can be measured as
        # the number of conflict grid cells with the same ID.
        # We could take the probability that a given conflict
        # spreads to remote cells as proportional to this size,
        # because there is more "coverage" of larger conflicts.
        #-------------------------------------------------------
        n_start = self.n_conflict_cells
        ID_count_max = np.max( ID_counts )
        
        for k in range(n_IDs):
            # Use U to exclude ocean & lake grid cells, etc.
            w2 = np.logical_and(self.S == 0, self.U != 0)
            ## w2 = (self.S == 0)   # (no-conflict cell locations)
            n2 = w2.sum()

            #----------------------------------------------            
            # Generate grid of Bernoulli random variables
            # using a grid of probabilities
            #----------------------------------------------
            if (n2 > 0):
                p = np.zeros(self.C2.shape, dtype='float32')
                p[ w2 ] = self.c_spread2 * ID_counts[k] / ID_count_max
                B  = np.random.binomial(1, p)  # entire grid
                ## B  = B.astype( 'uint8' )  ######
                w3 = (B == 1)   # (boolean array)
                n3 = w3.sum()

                self.S[ w3 ] = 1
                self.IDs[ w3 ]  = ID_vals[k]
                self.n_conflict_cells += n3

            #--------------------------------------
            # Alternate method that uses U and C2
            #--------------------------------------
            # not implemented
                                  
        #-----------------------------------------------------        
        # Compute probability that conflict in each cell with
        # conflict will spread to any other grid cell
        #-------------------------------------------------------
        # C2 measures how well-connected each grid cell is,
        # (via phone or internet access -- nonlocal) to other
        # grid cells (i.e. "the outside world"). This is
        # assumed proportional to the likelihood that conflict
        # will spread to them from distant grid cells.
        #-------------------------------------------------------      
        # w2 = np.invert(w1)   # (no-conflict cell locations)
        # C2_w2_max = np.max( self.C2[w2] ) 
        # p = np.zeros(self.C2.shape, dtype='float32')
        # p[ w2 ] = self.c_spread2 * self.C2[w2] / C2_w2_max

        #---------------------------------------------
        # Use Bernoulli r.v.s to determine spreading
        #---------------------------------------------
#         n_start = self.n_conflict_cells
#         
#         if (USE_LOOP):        
#             for k in range(n1):
#                 w2 = (self.S == 0)   # (no-conflict cell locations)
#                 C2_w2_max = np.max( self.C2[w2] ) 
#                 p = np.zeros(self.C2.shape, dtype='float32')
#                 if (C2_w2_max > 0):
#                     p[ w2 ] = self.c_spread2 * self.C2[w2] / C2_w2_max
# 
#                 #--------------------------------------
#                 # Alternate method that uses U and C2
#                 #--------------------------------------
#                 # not implemented
#                 
#                 # Generate grid of Bernoulli random variables
#                 B  = np.random.binomial(1, p)  # entire grid
#                 ## B  = B.astype( 'uint8' )  ######
#                 w3 = (B == 1)   # (boolean array)
#                 n3 = w3.sum()
# 
#                 self.S[ w3 ] = 1
#                 self.IDs[ w3 ]  = ID_vals[k]
#                 self.n_conflict_cells += n3
#         else:
#             pass
            #-------------------------------------
            # Spread conflict without a for loop
            # Much faster, and results look similar
            # See notes below re: overcounting.
            #-------------------------------------
#             B  = np.random.binomial(1, pn)  # (8 r.v.s)
#             ## B  = B.astype( 'uint8' )  ######
#             w2 = np.logical_and( (self.S[rn, cn] == 0), (B == 1) )
#             n2 = w2.sum()
#             #-----------------------------------------
#             # Spread conflict to some neighbor cells
#             #-----------------------------------------
#             # w2 is boolean array, but this is okay
#             #-----------------------------------------
#             self.S[ rn[w2], cn[w2] ]   = 1
#             self.IDs[ rn[w2], cn[w2] ] = ID_vals[w2]
#             #--------------------------------------------
#             # Without the for loop, several cells can 
#             # spread conflict to the same cell and this
#             # next line results in over-counting:
#             #    self.n_conflict_cells += n2
#             #--------------------------------------------
#             w_new = (self.S == 1)
#             n_new = (w_new.sum() - n1)
#             self.n_conflict_cells += n_new
                     
        if (self.REPORT):
            n_spread = (self.n_conflict_cells - n_start)
            print('Number of spread conflicts =', n_spread)
            print()
        
    #   spread_conflicts_nonlocal2()
    #---------------------------------------------------------------
    def finalize( self ):

        #-------------------------    
        # Close the output files
        #-------------------------
        self.out_unit.close()
        self.IDs_unit.close()
        
        if (self.REPORT):
            print()
            run_time = (time.time() - self.start_time)
            run_time = (run_time / 60.0)
            print('run_time  =', run_time, ' [minutes]')
            print('n_steps   =', self.n_steps)
            print('c_emerge  =', self.c_emerge, 'in (0,1)')
            print('c_spread  =', self.c_spread, 'in (0,1)')
            print('c_spread2 =', self.c_spread2, 'in (0,1)')
            ## print('p_geom   =', self.p_geom)
            print('p_resolve =', self.p_resolve, 'in (0,1)')
            print('time_lag  =', self.time_lag)
            #------------------------------------------
            # S has type 'uint8', so sum may not work
            #------------------------------------------
            ### n_conflict_cells = self.S.sum()
            w = (self.S == 1)
            n_conflict_cells = w.sum()
            #------------------------------------------
            print('n_conflict_cells =', n_conflict_cells)
            print('SELF CHECK...')
            print('n_conflict_cells =', self.n_conflict_cells)
            ## Exclude borders to count n_cells
            n_cells = (self.n_cols - 1) * (self.n_rows - 1)
            f_conflict_cells = (self.n_conflict_cells / n_cells)
            print('fraction_conflict_cells =', f_conflict_cells)
            print('Finished.')

    #   finalize()
    #---------------------------------------------------------------
    def run_model( self, cfg_file=None ):

        self.initialize( cfg_file )
        for k in range( self.n_steps ):
            self.update()
        self.finalize()
        
    #   run_model()
#-------------------------------------------------------------------


