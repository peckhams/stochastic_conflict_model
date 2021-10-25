
#  Copyright (c) 2021, Scott D. Peckham
#  May 2021.  Started on May 17 for Localized Conflict Modeling.
#             First, read GPW population count data as GeoTIFF.
#             Write-up and theoretical results in late May.
#-------------------------------------------------------------------

#  get_raster_cellsize()
#  get_raster_bounds()
#  bounds_disjoint()
#  read_geotiff()     # can also create RTG and RTI files
#  regrid_geotiff()

#-------------------------------------------------------------------

import gdal, osr  ## ogr
import glob, sys
import os, os.path
from . import rti_files
from . import rtg_files

#-------------------------------------------------------------------
def get_raster_cellsize( gdal_unit ):

    geotransform = gdal_unit.GetGeoTransform()
    # ulx  = geotransform[0]
    xres = geotransform[1]
    # xrtn = geotransform[2]
    #-----------------------
    # uly  = geotransform[3]
    # yrtn = geotransform[4]  # (not yres !!)
    yres = geotransform[5]  # (not yrtn !!)
    
    return (xres, abs(yres))

#   get_raster_cellsize()
#-------------------------------------------------------------------
def get_raster_bounds( ds, VERBOSE=False):

    #-------------------------------------------------------------
    # Note:  The bounds depend on the map projection and are not
    # necessarily a Geographic bounding box of lons and lats.    
    #-------------------------------------------------------------
    # See:
    # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
    # and search on "geotransform".  An example of gdal.SetGeoTransform
    # gives: [xmin, pixel_size, 0, ymax, 0, -pixel_size].
    # Also says args are:
    # [ulx, xDist, rtnX, uly, yDist, rtnY]
    # This is consistent with information below.
    #-------------------------------------------------------------    
    # ulx = upper left x  = xmin
    # uly = upper left y  = ymax
    # lrx = lower right x = xmax
    # lry = lower right y = ymin
    #-----------------------------

    #----------------------------------------------------------  
    # Notice the strange order or parameters here is CORRECT.
    # It is not:  ulx, xres, xskew, uly, yres, yskew
    #----------------------------------------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    if (VERBOSE):
        print('ulx, uly   =', ulx, uly)
        print('lrx, lry   =', lrx, lry)
        print('xres, yres = ', xres, yres)
        print('xskew, yskew =', xskew, yskew)
        print('----------------------------------')

    #########################################################
    # Bounding box reported by gdal.info does not match
    # what the GES DISC website is saying.  The result is
    # that gdal.Warp gives all nodata values in output. 
    #########################################################
    return [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]
 
    #########################################################
    # Bounding box reported by gdal.info does not match
    # what the GES DISC website is saying.  Reversing lats
    # and lons like this doesn't fix the problem.
    #########################################################    
    ## return [lry, ulx, uly, lrx]
    
#   get_raster_bounds()
#------------------------------------------------------------------- 
def bounds_disjoint( bounds1, bounds2, VERBOSE=False):
 
    #-----------------------------------------------------------
    # Note.  Assume both bounds are in same spatial reference
    #        system (SRS), e.g. Geographic lons and lats.
    #------------------------------------------------------------------
    # https://gamedev.stackexchange.com/questions/586/
    # what-is-the-fastest-way-to-work-out-2d-bounding-box-intersection
    #------------------------------------------------------------------    
    b1_xmin = bounds1[0]
    b1_xmax = bounds1[2]
    b2_xmin = bounds2[0]
    b2_xmax = bounds2[2]
#     x_overlap1 = (b1_xmin < b2_xmin) and (b2_xmin < b1_xmax)
#     x_overlap2 = (b2_xmin < b1_xmin) and (b1_xmin < b2_xmax)
#     x_overlap  = (x_overlap1 or x_overlap2)
    
    b1_ymin = bounds1[1]
    b1_ymax = bounds1[3]
    b2_ymin = bounds2[1]
    b2_ymax = bounds2[3]
#     y_overlap1 = (b1_ymin < b2_ymin) and (b2_ymin < b1_ymax) 
#     y_overlap2 = (b2_ymin < b1_ymin) and (b1_ymin < b2_ymax)
#     y_overlap  = (y_overlap1 or y_overlap2)
#     return not(x_overlap and y_overlap)

    disjoint = (b2_xmin > b1_xmax) or (b2_xmax < b1_xmin) or \
               (b2_ymax < b1_ymin) or (b2_ymin > b1_ymax)

    return disjoint
    
#   bounds_disjoint()
#-------------------------------------------------------------------
def read_geotiff(in_file=None, REPORT=True,
                 MAKE_RTG=False, rtg_file=None ):
 
    # Bounds = [ minlon, minlat, maxlon, maxlat ]'

    if (in_file is None):
        in_dir  = '/Users/peckhams/Conflict/Data/GPW-v4/'  
        in_file = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
        in_file = in_dir + in_file

    if (rtg_file is None):
        in_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
        rtg_file = 'GPW-v4_global_pop_count.rtg'
        rtg_file = in_dir + rtg_file
          
    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    print('Reading grid from GeoTIFF file...')
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_xres_deg = dx
    in_yres_deg = dy
    in_xres_sec = (in_xres_deg * 3600.0)
    in_yres_sec = (in_yres_deg * 3600.0)
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    #------------------------------
    # Get min & max of input grid
    #------------------------------
    in_grid  = in_unit.ReadAsArray()
    in_dtype = in_grid.dtype
    in_gmin  = in_grid.min()
    in_gmax  = in_grid.max()
     
    #----------------       
    # Close in_file
    #----------------
    in_unit = None   # Close in_file

    #-----------------------------
    # Option to save as RTG file
    #-----------------------------
    if (MAKE_RTG):
        rtg_path = in_dir + rtg_file   #####
        #-------------------------------
        # Option to create an RTI file
        #-------------------------------
        rti = rti_files.make_info(grid_file=rtg_path,
                  ncols=in_ncols, nrows=in_nrows,
                  xres=in_xres_sec, yres=in_yres_sec,
                  #--------------------------------------
                  data_source='SEDAC',
                  data_type='FLOAT',
                  byte_order=rti_files.get_rti_byte_order(),
                  pixel_geom=0,
                  zres=0.001, z_units='unknown',
                  y_south_edge=in_bounds[1],
                  y_north_edge=in_bounds[3],
                  x_west_edge=in_bounds[0],
                  x_east_edge=in_bounds[2],
                  box_units='DEGREES')
        rti_files.write_info( rtg_path, rti )    
        rtg = rtg_files.rtg_file() 
        OK  = rtg.open_new_file( rtg_path, rti )
        if not(OK):
            print('ERROR during open_new_file().')
            return       
        rtg.write_grid( in_grid, VERBOSE=True )
        rtg.close_file()
         
    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('GeoTIFF grid info:')
        print('   ' + in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        print('   xres   =', in_xres_sec, ' [arcsecs]' )
        print('   yres   =', in_yres_sec, ' [arcsecs]' )
        print('   dtype  =', in_dtype )
        print('   bounds =', in_bounds )
        print('   gmin   =', in_gmin )
        print('   gmax   =', in_gmax ) 
        print('Finished.')
        print()

    #------------------------------
    # Return the grid as an array
    #------------------------------
    return in_grid

#   read_geotiff()
#-------------------------------------------------------------------
def regrid_geotiff(in_file=None, out_file=None, 
                   out_bounds=None,
                   out_xres_sec=None, out_yres_sec=None,
                   RESAMPLE_ALGO='bilinear', REPORT=True):
                   ##### in_nodata=None, out_nodata=None):
                   
    #--------------------------------------------------------
    # Notes:  Read grid from GeoTIFF file, optionally clip
    #         and resample it, then save result as GeoTIFF.
    #--------------------------------------------------------

    #-----------------------------------   
    # Specify the resampling algorithm
    #-----------------------------------
    algo_dict = {
    'nearest'     : gdal.GRA_NearestNeighbour,
    'bilinear'    : gdal.GRA_Bilinear,
    'cubic'       : gdal.GRA_Cubic,
    'cubicspline' : gdal.GRA_CubicSpline,
    'lanczos'     : gdal.GRA_Lanczos,
    'average'     : gdal.GRA_Average,
    'min'         : gdal.GRA_Min,
    'max'         : gdal.GRA_Max,
    'mode'        : gdal.GRA_Mode,
    'med'         : gdal.GRA_Med }
    
    resample_algo = algo_dict[ RESAMPLE_ALGO ]
    
    #---------------------------------------------------------------
    # Note:  DEM_bounds = [dem_xmin, dem_ymin, dem_xmax, dem_ymax]
    #        Give xres, yres in *arcseconds* for Geographic.
    #        They will then be converted to decimal degrees.
    #        gdal.Warp() clips to a bounding box, and can also
    #        resample to a different resolution.
    #        gdal.Translate() is faster for simple clipping.
    #---------------------------------------------------------------
    if (in_file == None):
        #--------------------------------
        # Use Horn of Africa as a test.
        #--------------------------------
        test_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
        in_file    = test_dir + 'gpw_v4_population_count_rev11_2020_30_sec.tif'
        out_file   = test_dir + 'Horn_of_Africa_pop_count.tif'
        out_xres_sec = None    # will default to in_xres_sec
        out_yres_sec = None    # will default to in_yres_sec
        # Bounds = [ minlon, minlat, maxlon, maxlat ]
        out_bounds = [ 25.0, -5.0, 55.0, 25.0]
  
    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_xres_deg = dx
    in_yres_deg = dy
    in_xres_sec = (in_xres_deg * 3600.0)
    in_yres_sec = (in_yres_deg * 3600.0)
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    #------------------------------
    # Get min & max of input grid
    #------------------------------
    in_grid  = in_unit.ReadAsArray()
    in_dtype = in_grid.dtype
    in_gmin  = in_grid.min()
    in_gmax  = in_grid.max()
    
    #-------------------------------------------------------------
    # If a spatial resolution has not been specified for output,
    # then assume it is the same as the resolution of the input.
    #-------------------------------------------------------------
    if (out_xres_sec is not None):
        out_xres_deg = (out_xres_sec / 3600.0)  # [arcsec] -> [degrees]
    else:
        out_xres_deg = None    # (will default to in_xres_deg)
        ## out_xres_sec = in_xres_sec
        ## out_xres_deg = (out_xres_sec / 3600.0)  # [arcsec] -> [degrees]       
    #----------------------------------------------------------------------
    if (out_yres_sec is not None):
        out_yres_deg = (out_yres_sec / 3600.0)  # [arcsec] -> [degrees]
    else:
        out_yres_deg = None  # (will default to in_yres_deg)
        ## out_yres_sec = in_yres_sec
        ## out_yres_deg = (out_yres_sec / 3600.0)  # [arcsec] -> [degrees]     
                
    #------------------------------------------- 
    # Are out_bounds disjoint from in_bounds ?
    #-------------------------------------------
    # Specifying out_bounds is not required.
    #-------------------------------------------
    if (out_bounds is not None):   
        DISJOINT = bounds_disjoint( in_bounds, out_bounds, VERBOSE=False)
        if (DISJOINT):
            # print('ERROR:  Output bounds do not overlap input bounds.')
            print('ERROR:  Input & output bounding boxes are disjoint.')
            print( '       New grid would contain only nodata.')
            print('  in_file  =', in_file )
            print('  out_file =', out_file )
            in_unit  = None   # Close in_file
            return

    #------------------------------------------------  
    # Resample & clip and write new grid to GeoTIFF
    #------------------------------------------------
    out_unit = gdal.Warp( out_file, in_unit,
        format = 'GTiff',  # (output format string)
        outputBounds=out_bounds,
        xRes=out_xres_deg, yRes=out_yres_deg,
        # srcNodata = in_nodata,      ########  FUTURE
        # dstNodata = out_nodata,     ########  FUTURE
        resampleAlg = resample_algo )

    #-----------------------
    # Get info on new grid
    #-----------------------
    out_ncols  = out_unit.RasterXSize
    out_nrows  = out_unit.RasterYSize
    out_bounds = get_raster_bounds( out_unit )   ######

    #-------------------------------
    # Get min & max of output grid
    #-------------------------------
    out_grid  = out_unit.ReadAsArray()
    out_dtype = out_grid.dtype
    out_gmin  = out_grid.min()
    out_gmax  = out_grid.max()
        
    #----------------------------------        
    # Close both in_file and out_file
    #----------------------------------
    in_unit  = None   # Close in_file
    out_unit = None   # Close out_file

    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('Input grid file:')
        print('   ' + in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        print('   xres   =', in_xres_sec, ' [arcsecs]' )
        print('   yres   =', in_yres_sec, ' [arcsecs]' )
        print('   bounds =', in_bounds )
        print('   dtype  =', in_dtype )
        print('   gmin   =', in_gmin )
        print('   gmax   =', in_gmax )
        print()
        print('Output grid file:')
        print('   ' + out_file )
        print('   ncols  =', out_ncols )
        print('   nrows  =', out_nrows )
        print('   xres   =', out_xres_sec, ' [arcsecs]' )
        print('   yres   =', out_yres_sec, ' [arcsecs]' ) 
        print('   bounds =', out_bounds )
        print('   dtype  =', out_dtype )
        print('   gmin   =', out_gmin )
        print('   gmax   =', out_gmax )      
        print('Finished regridding.')
        print()

    #--------------------------------------------------------  
    # This shows some of the other keywords to gdal.Warp.
    #--------------------------------------------------------      
    # WarpOptions(options=[], format=None, outputBounds=None,
    # outputBoundsSRS=None, xRes=None, yRes=None, targetAlignedPixels=False,
    # width=0, height=0, srcSRS=None, dstSRS=None, srcAlpha=False,
    # dstAlpha=False, warpOptions=None, errorThreshold=None,
    # warpMemoryLimit=None, creationOptions=None, outputType=GDT_Unknown,
    # workingType=GDT_Unknown, resampleAlg=None, srcNodata=None,
    # dstNodata=None, multithread=False, tps=False, rpc=False,
    # geoloc=False, polynomialOrder=None, transformerOptions=None,
    # cutlineDSName=None, cutlineLayer=None, cutlineWhere=None,
    # cutlineSQL=None, cutlineBlend=None, cropToCutline=False,
    # copyMetadata=True, metadataConflictValue=None,
    # setColorInterpretation=False, callback=None, callback_data=None)
    # 
    # Create a WarpOptions() object that can be passed to gdal.Warp()
    # Keyword arguments are : options --- can be be an array of strings,
    # a string or let empty and filled from other keywords.

#    regrid_geotiff()
#-------------------------------------------------------------------



