#   
#  Copyright (c) 2020, Scott D. Peckham
#
#  Note: This file contains a set of functions for visualizing the
#        contents of output files in netCDF format
#        (e.g. TopoFlow or Stochastic Conflict Model)
#
#  May 2020.  Moved all routines from Jupyter notebook called
#             TopoFlow_Visualization.ipynb to here.
#             Tested all of them in a new Jupyter notebook called
#             TopoFlow_Visualization2.ipynb.
#
#--------------------------------------------------------------------
#
#  Define some stretch functions for 2D color images:
#  normalize_grid()
#  histogram_equalize()
#  power_stretch0()
#  power_stretch1()
#  power_stretch2()
#  power_stretch3()
#  log_stretch()
#  linear_stretch()
#  stretch_grid()
#
#  Define functions to show grids as color images:
#  read_grid_from_nc_file()
#  read_and_show_rtg()
#  show_grid_as_image()
#  save_grid_stack_as_images()
#  save_rts_as_images()
#
#  Create movies from set of images:
#     (works for grid images, profile images, etc.)
#  create_movie_from_images()
#
#  plot_data()
#  create_media_files()
#  delete_png_files()

#--------------------------------------------------------------------
# import os.path
# import shutil

import glob, os, time
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import imageio
import imageio_ffmpeg as imff  # for adding opacity

from conflict.utils import ncgs_files
from conflict.utils import rtg_files
from conflict.utils import rts_files

#--------------------------------------------------------------------
def normalize_grid( grid ): 

    gmin = grid.min()
    gmax = grid.max()

    if (gmin != gmax):
        norm = (grid - gmin) / (gmax - gmin)
    else:
        # Avoid divide by zero
        norm = np.zeros( grid.shape, dtype=grid.dtype )
    return norm

#   normalize_grid()
#--------------------------------------------------------------------
def histogram_equalize( grid, PLOT_NCS=False):

    #  https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
    (hist, bin_edges) = np.histogram( grid, bins=256)
    # hmin = hist.min()
    # hmax = hist.max()

    cs  = hist.cumsum()
    ncs = (cs - cs.min()) / (cs.max() - cs.min())
    ncs.astype('uint8');

    if (PLOT_NCS):
        plt.plot( ncs )

    flat = grid.flatten()
    if (flat.max() != flat.min()):
        flat2 = np.uint8( 255 * (flat - flat.min()) / (flat.max() - flat.min()) )
        grid2 = ncs[ flat2 ].reshape( grid.shape )
    else:
        flat2 = np.zeros( flat.size, dtype='uint8' )
        grid2 = ncs[ flat2 ].reshape( grid.shape )

    return grid2

#   histogram_equalize()
#--------------------------------------------------------------------
def power_stretch0( grid, p ):

    norm = normalize_grid( grid )
    
    return norm**p
    
#   power_stretch0()
#--------------------------------------------------------------------
def power_stretch1( grid, p ):
    return grid**p
    
#   power_stretch1()
#--------------------------------------------------------------------
def power_stretch2( grid, a=1000, b=0.5):

    # Note: Try a=1000 and b=0.5
    norm = normalize_grid( grid )
    return (1 - (1 + a * norm)**(-b))
    
#   power_stretch2()
#--------------------------------------------------------------------
def power_stretch3( grid, a=1, b=2):

    # Note:  Try a=1, b=2 (shape of a quarter circle)
    norm = normalize_grid( grid )
    return (1 - (1 - norm**a)**b)**(1/b)
    
#   power_stretch3()
#--------------------------------------------------------------------
def log_stretch( grid, a=1 ):
    return np.log( (a * grid) + 1 )
    
#   log_stretch()
#--------------------------------------------------------------------
def linear_stretch( grid ):

    norm = normalize_grid( grid )
    return norm
   
#   linear_stretch()
#--------------------------------------------------------------------
def stretch_grid( grid, stretch, a=1, b=2, p=0.5 ):

    name = stretch
    if   (name == 'hist_equal'):
        grid2 = histogram_equalize( grid, PLOT_NCS=False)    
    elif (name == 'linear'):
        grid2 = linear_stretch(grid)
    elif (name == 'log'):
        grid2 = log_stretch( grid, a=a )
    elif (name == 'power'):
        grid2 = power_stretch0( grid, p=p )
    elif (name == 'power1'): 
        # Try:  p = 0.3   
        grid2 = power_stretch1( grid, p)
    elif (name == 'power2'):
        # Try:  a=1000, b=0.5.
        grid2 = power_stretch2( grid, a=a, b=b )
    elif (name == 'power3'):        
        # Try:  a=1, b=2.
        grid2 = power_stretch3( grid, a=a, b=b)
    else:
        print('### SORRY, Unknown stretch =', name)
        return grid

    return grid2
 
#   stretch_grid()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def read_grid_from_nc_file( nc_file, time_index=1, REPORT=True ):

    # Typical 2D nc files
    # nc_file = case_prefix + '_2D-Q.nc'
    # nc_file = case_prefix + '_2D-d-flood.nc'

    ## if ('_2D' not in nc_file):
    if ('2D' not in nc_file):
        print('ERROR: This function is only for "2D" netCDF files.')
        print('       Filename must contain "2D" for grid stacks.') 
        return
            
    ncgs = ncgs_files.ncgs_file()
    ncgs.open_file( nc_file )
    var_name_list = ncgs.get_var_names()
    if (REPORT):
        print('var_names in netCDF file =' )
        print( var_name_list )

    #----------------------------         
    # Determine valid var_index
    #-----------------------------------------
    # Old: 0=time, 1=X, 2=Y, 3=V
    # New: 0=time, 1=datetime, 2=X, 3=Y, 4=V
    #-----------------------------------------
    var_index = 1
    other_vars = ['time','datetime','X','Y','Z']
    while (True):
        var_name = var_name_list[ var_index ]
        if (var_name not in other_vars):
            break
        var_index += 1    
    ### var_index = 3   # 0=time, 1=X, 2=Y, 3=V  ###############
    ### var_name  = var_name_list[ var_index ]
    long_name = ncgs.get_var_long_name( var_name )
    var_units = ncgs.get_var_units( var_name )
    n_grids   = ncgs.ncgs_unit.variables[ var_name ].n_grids

    if (REPORT):
        print('long_name =', long_name)
        print('var_name  =', var_name)
        print('var_units =', var_units)
        print('n_grids   =', n_grids)

    #--------------------------------------------
    # Use these to set "extent" in plt.imshow()
    #--------------------------------------------
    minlon = ncgs.ncgs_unit.variables['X'].geospatial_lon_min
    maxlon = ncgs.ncgs_unit.variables['X'].geospatial_lon_max
    minlat = ncgs.ncgs_unit.variables['Y'].geospatial_lat_min
    maxlat = ncgs.ncgs_unit.variables['Y'].geospatial_lat_max
    extent = [minlon, maxlon, minlat, maxlat]
    
    #----------------------------------------------
    # Read grid from nc_file for given time_index
    #----------------------------------------------
    grid = ncgs.get_grid( var_name, time_index )
    
    if (REPORT):
        print( 'extent = ')
        print( extent )
        print( 'grid shape =', grid.shape )
        print( 'min(grid)  =', grid.min() )
        print( 'max(grid)  =', grid.max() )

    ncgs.close_file()
    return (grid, long_name, extent)
    
#   read_grid_from_nc_file()
#--------------------------------------------------------------------
def read_and_show_rtg( rtg_filename, long_name, VERBOSE=True,
                       cmap='jet', BLACK_ZERO=False,
                       stretch='hist_equal',
                       a=1, b=2, p=0.5, im_file=None,
                       xsize=8, ysize=8, dpi=None ):
    
    rtg = rtg_files.rtg_file()
    OK  = rtg.open_file( rtg_filename )
    if not(OK):
        print('Sorry, Could not open RTG file:')
        print( rtg_filename )
        return
    
    grid   = rtg.read_grid( VERBOSE=VERBOSE )
    extent = rtg.get_bounds()
    rtg.close_file()

    if (VERBOSE):
        print('Byte swap needed =', rtg.byte_swap_needed())
        print('Reading grid from RTG file...')
        print('extent =', extent)
        print('min(grid), max(grid) =', grid.min(), grid.max())
        print('Finished.')
        print()

    show_grid_as_image( grid, long_name, extent=extent, cmap=cmap,
                        BLACK_ZERO=BLACK_ZERO, stretch=stretch,
                        a=a, b=b, p=p, im_file=im_file,
                        xsize=xsize, ysize=ysize, dpi=dpi)
                              
#   read_and_show_rtg()
#--------------------------------------------------------------------
def show_grid_as_image( grid, long_name, extent=None,
                        cmap='rainbow', BLACK_ZERO=False,
                        LAND_SEA_BACKDROP=False,
                        LAND_SEA_RED_BACKDROP=False,                        
                        stretch='power3',
                        a=1, b=2, p=0.5,
                        NO_SHOW=False, im_file=None,
                        xsize=8, ysize=8, dpi=None): 

    # Note:  extent = [minlon, maxlon, minlat, maxlat]
    
    #-------------------------
    # Other color map names
    #--------------------------------------------
    # hsv, jet, gist_rainbow (reverse rainbow),
    # gist_ncar, gist_stern
    #--------------------------------------------    

    #--------------------------
    # Other stretch functions
    #--------------------------
    grid2 = stretch_grid( grid, stretch, a=a, b=b, p=p )
#     if (stretch == 'power_stretch3'):
#         grid2 = power_stretch3( grid, a=0.5 )
#     elif (stretch == 'power_stretch1a'):   
#         grid2 = power_stretch1( grid, 0.5)
#     elif (stretch == 'power_stretch1b'):
#         grid2 = power_stretch1( grid, 0.2)
#     elif (stretch == 'power_stretch2'):
#         grid2 = power_stretch2( grid )
#     elif (stretch == 'log_stretch'):
#         grid2 = log_stretch( grid )
#     elif (stretch == 'hist_equal'):
#         grid2 = histogram_equalize( grid, PLOT_NCS=True)
#     else:
#         print('SORRY, Unknown stretch =', stretch)
#         return 

    #---------------------------------------
    # Modify the colormap (0 = black) ?
    # cmap is name of colormap, a string
    #--------------------------------------------------------
    # cmap arg to imshow can be name (as str) or cmap array
    # 4th entry is opacity, or alpha channel (I think)
    #--------------------------------------------------------
    # See: "Creating listed colormaps" section at:
    # https://matplotlib.org/3.1.0/tutorials/colors/
    #         colormap-manipulation.html
    #--------------------------------------------------------
    # "Land green" = #c6e5bc = (198, 229, 188)
    # "Sea blue"   = #aad3df = (170, 211, 223)
    #--------------------------------------------------------
    if (BLACK_ZERO):
        n_colors = 256
        color_map  = cm.get_cmap( cmap, n_colors )
        new_colors = color_map( np.linspace(0, 1, n_colors) )
        black = np.array([0.0, 0.0, 0.0, 1.0])
        new_colors[0,:] = black
        new_cmap = ListedColormap( new_colors )
    elif (LAND_SEA_BACKDROP):
        n_colors = 256
        color_map  = cm.get_cmap( cmap, n_colors )
        new_colors = color_map( np.linspace(0, 1, n_colors) )
        land_green = np.array([198, 229, 188, 256]) / 256.0
        sea_blue   = np.array([170, 211, 223, 256]) / 256.0
        new_colors[0,:]   = land_green
        new_colors[255,:] = sea_blue
        new_cmap = ListedColormap( new_colors )
    elif (LAND_SEA_RED_BACKDROP):
        n_colors = 3
        color_map  = cm.get_cmap( cmap, n_colors )
        new_colors = color_map( np.linspace(0, 1, n_colors) )
        land_green = np.array([198, 229, 188, 256]) / 256.0
        sea_blue   = np.array([170, 211, 223, 256]) / 256.0
        red        = np.array([256,   0,   0, 256]) / 256.0
        new_colors[0,:] = land_green
        new_colors[1,:] = red
        new_colors[2,:] = sea_blue
        new_cmap = ListedColormap( new_colors )
    else:
        new_cmap = cmap
    
    #----------------------------
    # Set up and show the image
    #----------------------------
    # figure = plt.figure(1, figsize=(xsize, ysize))
    fig, ax = plt.subplots( figsize=(xsize, ysize), dpi=dpi)
    im_title = long_name.replace('_', ' ').title()
    #----------------------------------------------------
    if (dpi is not None):
        plt.rc('axes',  titlesize=8) # fontsize of the title
        plt.rc('axes',  labelsize=6)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=6)  # fontsize of the x tick labels
        plt.rc('ytick', labelsize=6)  # fontsize of the y tick labels
        # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    #----------------------------------------------------
    # Default label size is 10 points    
    ax.set_title( im_title)
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    #----------------------------------------------------
    # plt.setp(ax.get_xticklabels(), fontsize=8)
    # plt.setp(ax.get_yticklabels(), fontsize=8)
    #----------------------------------------------------
#     ax.set_title( im_title, fontsize=16 )
#     ax.set_xlabel('Longitude [deg]', fontsize=14)
#     ax.set_ylabel('Latitude [deg]',  fontsize=14)
#     ax.tick_params(axis='x', labelsize=14)
#     ax.tick_params(axis='y', labelsize=14)
    #---------------------------------------
    gmin = grid2.min()
    gmax = grid2.max()
    im = ax.imshow(grid2, interpolation='nearest', cmap=new_cmap,
                   vmin=gmin, vmax=gmax, extent=extent)

    #--------------------------------------------------------        
    # NOTE!  Must save before "showing" or get blank image.
    #        File format is inferred from extension.
    #        e.g. TMP_Image.png, TMP_Image.jpg.
    #--------------------------------------------------------
    if (im_file is not None):  
        plt.savefig( im_file )
    else:
        plt.show()   # Ignore NO_SHOW arg for now.   
    #-----------------------------------------------
#     if (im_file is not None):  
#         plt.savefig( im_file )
#     if not(NO_SHOW):
#         plt.show()
 
    plt.close()
        
#   Information on matplotlib color maps
#   https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
# 
#   Information on matplotlib.pyplot.imshow
#   https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.imshow.html
# 
#   Information on matplotlib.pyplot.savefig
#   https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.savefig.html
# 
#   plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
#               orientation='portrait', papertype=None, format=None,
#               transparent=False, bbox_inches=None, pad_inches=0.1,
#               frameon=None, metadata=None)
  
#   show_grid_as_image()
#--------------------------------------------------------------------
def save_grid_stack_as_images( nc_file, png_dir, extent=None,
              stretch='power3', a=1, b=2, p=0.5,
              cmap='jet', REPORT=True,
              ## cmap='rainbow', REPORT=True,
              BLACK_ZERO=False,
              LAND_SEA_BACKDROP=False,
              LAND_SEA_RED_BACKDROP=False,
              xsize=6, ysize=6, dpi=None ):
              ### xsize=6, ysize=6, dpi=192 ):

    #----------------------------------------------------
    # Note: dpi=None gives good default fontsizes, etc.
    #----------------------------------------------------

    ## if ('_2D' not in nc_file):
    if ('2D' not in nc_file):
        print('ERROR: This function is only for "2D" netCDF files.')
        print('       Filename must contain "2D" for grid stacks.') 
        return

    ncgs = ncgs_files.ncgs_file()        
    ncgs.open_file( nc_file )
    var_name_list = ncgs.get_var_names( no_dim_vars=True )  ####
    var_index = 0   # (dim vars are now excluded)
    var_name  = var_name_list[ var_index ]
    long_name = ncgs.get_var_long_name( var_name )
    n_grids = ncgs.ncgs_unit.variables[var_name].n_grids

    #-----------------------------------------------
    # MINT netCDF conventions:
    #    bounds = [minlon, minlat, maxlon, maxlat]
    # For matplotlib.pyplot.imshow():
    #    extent = [minlon, maxlon, minlat, maxlat]
    #-----------------------------------------------
    if (extent is None):
        bounds = ncgs.ncgs_unit.geospatial_bounds 
        extent = [ bounds[0], bounds[2], bounds[1], bounds[3]]

    im_title = long_name.replace('_', ' ').title()
    im_file_prefix = 'TF_Grid_Movie_Frame_'
    time_pad_map = {1:'0000', 2:'000', 3:'00', 4:'0', 5:''}
    cmap = 'rainbow'

    if (REPORT):
        print('Creating images from grid stack in nc_file:')
        print('  ' + nc_file )
        print('  ' + 'var name  =', var_name)
        print('  ' + 'long name =', long_name)
        print('  ' + 'n_grids   =', n_grids)
        print('This may take a few minutes.')
        print('Working...')
   
    time_index = 0
    while (True):
        # print('time index =', time_index )
        try:
            grid = ncgs.get_grid( var_name, time_index )
        except:
            break
        time_index += 1

        #----------------------------------------    
        # Build a filename for this image/frame
        #----------------------------------------
        tstr = str(time_index)
        pad = time_pad_map[ len(tstr) ]
        time_str = (pad + tstr)
        im_file = im_file_prefix + time_str + '.png' 
        im_file = (png_dir + '/' + im_file)
                
        show_grid_as_image( grid, long_name, cmap=cmap,
                            stretch=stretch, a=a, b=b, p=p, 
                            extent=extent,
                            BLACK_ZERO=BLACK_ZERO,
                            LAND_SEA_BACKDROP=LAND_SEA_BACKDROP,
                            LAND_SEA_RED_BACKDROP=LAND_SEA_RED_BACKDROP,
                            NO_SHOW=True, im_file=im_file,
                            xsize=xsize, ysize=ysize, dpi=dpi )
                            
    ncgs.close_file()
    tstr = str(time_index)
    print('Finished saving ' + tstr + ' images to PNG files.')

#   save_grid_stack_as_images()
#--------------------------------------------------------------------
def save_rts_as_images( rts_file, png_dir, extent=None,
                        long_name='River Discharge',
                        stretch='power3', a=1, b=2, p=0.5,
                        cmap='rainbow', BLACK_ZERO=False,
                        REPORT=True,
                        xsize=6, ysize=6, dpi=None):
                        ### xsize=6, ysize=6, dpi=192):

    #----------------------------------------------------
    # Note: dpi=None gives good default fontsizes, etc.
    #----------------------------------------------------

    if ('.rts' not in rts_file):
        print('ERROR: This function is only for RTS files.') 
        return

    rts = rts_files.rts_file()
    OK  = rts.open_file( rts_file )
    if not(OK):
        print('Could not open RTS file.')
        return
    n_grids = rts.number_of_grids()
    print('Byte swap needed =', rts.byte_swap_needed())

    # For matplotlib.pyplot.imshow():
    # extent = [minlon, maxlon, minlat, maxlat]    
    if (extent is None):
        extent = rts.get_bounds()

    im_title = long_name.replace('_', ' ').title()
    im_file_prefix = 'TF_RTS_Movie_Frame_'
    time_pad_map = {1:'0000', 2:'000', 3:'00', 4:'0', 5:''}

    if (REPORT):
        print('Creating images from grid stack in rts_file:')
        print('  ' + rts_file )
        print('  ' + 'long name =', long_name)
        print('  ' + 'n_grids   =', n_grids)
        print('  ' + 'extent    =', extent)
        print('This may take a few minutes.')
        print('Working...')
        
    time_index = 0
    rts_min = 1e12
    rts_max = 1e-12

    while (True):
        # print('time index =', time_index )
        try:
            grid = rts.read_grid( time_index )   # alias to get_grid()
            gmin = grid.min()
            gmax = grid.max()
            rts_min = min( rts_min, gmin )
            rts_max = max( rts_max, gmax )
        except:
            break
        time_index += 1

        #----------------------------------------    
        # Build a filename for this image/frame
        #----------------------------------------
        tstr = str(time_index)
        pad = time_pad_map[ len(tstr) ]
        time_str = (pad + tstr)
        im_file = im_file_prefix + time_str + '.png' 
        im_file = (png_dir + '/' + im_file)
                
        show_grid_as_image( grid, long_name, cmap=cmap,
                            stretch=stretch, a=a, b=b, p=p,
                            BLACK_ZERO=BLACK_ZERO, extent=extent,
                            NO_SHOW=True, im_file=im_file,
                            xsize=xsize, ysize=ysize, dpi=dpi )

    rts.close_file()
    print('min(rts), max(rts) =', rts_min, rts_max)
    tstr = str(time_index)
    print('Finished saving ' + tstr + ' images to PNG files.')  
    print()

#   save_rts_as_images()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def create_movie_from_images( mp4_file, png_dir, fps=10, REPORT=True):

    #----------------------------------------
    # png_dir  = directory with PNG files
    # mp4_file = case_prefix + '_Movie.mp4'
    # fps      = frames per second
    #----------------------------------------
    im_file_list = sorted( glob.glob( png_dir + '/*.png' ) )
    n_frames = len( im_file_list )

    if (REPORT):
        print('Creating movie from', n_frames, 'PNG files.')
        ## print('This may take a few minutes.')
        print('Working...')

    #---------------------------------------------------------------        
    # Note:  The default codec for imageio is "libx264", which
    #        is widely used and supports mp4; H.264 codec.
    #        You can request  "mpeg4" also, and it works.
    #        The codec "mjpeg" is said to be good for JPEG images.
    #        If you use Get Info on an MP4, the More Info section
    #        give codecs for these as:  "H.264" or "MPEG-4 Video".
    #        If I copy an MP4 to: topoflow36/movies on GitHub,
    #        I can't get them to play in Jupyter notebook or lab.
    #        See the notebook:  MP4_Video_Test.ipynb for more info.
    # https://imageio.readthedocs.io/en/stable/format_ffmpeg.html
    # https://imageio.readthedocs.io/en/stable/reference/_backends/
    #         imageio.plugins.ffmpeg.html
    #----------------------------------------------------------------
    #  quality ranges from 0 to 10.  Default is 5 in ffmpeg plugin.
    #    But QuickTime can't play if quality = 10 & default codec.
    #  pixelformat = 'yuvj444p' records all 3 channels at full res.
    #                same as 'yuv444p' ???
    #                ffmpeg plugin default is: 'yuv420p'
    #----------------------------------------------------------------
    #writer = imageio.get_writer( mp4_file, fps=fps, quality=9.5 ) # works 1.2-1.9MB (IDs)
    # writer = imageio.get_writer( mp4_file, fps=fps, quality=9)  # works 1.5MB (IDs)
    # writer = imageio.get_writer( mp4_file, fps=fps, quality=8)  # works 1.4MB (IDs)
    writer = imageio.get_writer( mp4_file, fps=fps )  # (default quality=5)
    # writer = imageio.get_writer( mp4_file, fps=fps, codec='mpeg4' )

    #------------------------------
    # QuickTime can't play these.
    #------------------------------
    # writer = imageio.get_writer( mp4_file, fps=fps, quality=10 )  # 1.3MB (IDs)
    # writer = imageio.get_writer( mp4_file, fps=fps, quality=9.9 ) # 1.3MB (IDs)
    # writer = imageio.get_writer( mp4_file, fps=fps, quality=9.6 ) # 1.6MB (IDs)
    # writer = imageio.get_writer( mp4_file, fps=fps,
    #                  quality=10, pixelformat='yuvj444p' )
    # writer = imageio.get_writer( mp4_file, fps=fps,
    #                  quality=9, pixelformat='yuvj444p' )
    # writer = imageio.get_writer( mp4_file, fps=fps,
    #                  quality=8, pixelformat='yuvj444p' )
    # writer = imageio.get_writer( mp4_file, fps=fps,
    #                  quality=8, pixelformat='yuv444p' )

    for im_file in im_file_list:
        writer.append_data(imageio.imread( im_file ))
    writer.close()
   
    if (REPORT): 
        print('Finished creating movie, MP4 format.')
        print('  ' + mp4_file)
        print()

#   create_movie_from_images()
#--------------------------------------------------------------------
def plot_data( x, y, y2=None, xmin=None, xmax=None, ymin=None, ymax=None,
               x_name='x', x_units='', marker=',', title=None,
               y_name='y', y_units='',
               x_size=8,   y_size=4):

    figure = plt.figure(1, figsize=(x_size, y_size))
    # fig, ax = plt.subplots( figsize=(x_size, y_size))

    # Set the plot point marker
    # https://matplotlib.org/3.1.1/api/markers_api.html
    # marker = ','  # pixel
    # marker = '.'  # point (small circle)
    # marker = 'o'  # circle
    # marker = '+'
    # marker = 'x'

    #if (ymin is None):
    #    ymin = y.min()
    #if (ymax is None):
    #    ymax = y.max()
    #if (ymax - ymin < 0.1):
    #    ymin = ymin - 0.5
    #    ymax = ymin + 0.5

    # x_name2 = x_name.replace('_', ' ').title()
    # y_name2 = y_name.replace('_', ' ').title()
        
    plt.plot( x, y, marker=marker)
    if (y2 is not None):
        plt.plot(x, y2, marker=marker)

    plt.xlabel( x_name + ' [' + x_units + ']' )
    plt.ylabel( y_name + ' [' + y_units + ']' )
    if (title is not None):
        plt.title( title )

    plt.ylim( ymin, ymax )
    plt.xlim( xmin, xmax )
    #-------------------------------------
    # This may be necessary depending on
    # the data type of ymin, ymax
    #-------------------------------------
    ## plt.ylim( np.array([ymin, ymax]) )
    ## plt.xlim( np.array([xmin, xmax]) )
    plt.show()

#   plot_data()
#----------------------------------------------------------------------------   
def create_media_files( output_dir=None, media_dir=None,
                        WEBM=True, opacity=1.0,
                        xsize2D=4, ysize2D=4,
                        movie_fps=10, dpi=None):
                        ## movie_fps=10, dpi=192):

    #----------------------------------------------------
    # Note: dpi=None gives good default fontsizes, etc.
    #----------------------------------------------------
                                    
    #------------------------------------------------
    # Write a separate function to create movies
    # from the rainfall grid stacks (RTS format) ??
    #--------------------------------------------------
    # NOTE!  Dojo doesn't allow media directory to be
    #        in output directory.  Must be siblings
    #        because they are mounted separately.
    #        Still need output_dir; added media_dir.
    #--------------------------------------------------
    if (output_dir is None):
        print('SORRY, output_dir argument is required.')
        print()
        return
    
    if (media_dir is None):
        print('SORRY, media_dir argument is required.')
        print()
        return

    #-------------------------------------------
    # Get & save the current working directory
    # likely, the package directory
    #-------------------------------------------
    start_dir = os.getcwd()
                   
    #------------------------------------
    # Setup required output directories
    #------------------------------------
    temp_png_dir = media_dir + 'temp_png/'
    if not(os.path.exists( temp_png_dir )):
        os.mkdir( temp_png_dir )
    #-----------------------------------------
    movie_dir = media_dir + 'movies/'
    if not(os.path.exists( movie_dir )):
        os.mkdir( movie_dir )
    #-----------------------------------------
    # image_dir = media_dir + 'images/'
    # if not(os.path.exists( image_dir )):
    #     os.mkdir( image_dir )
  
    #----------------------------------------------
    # Create set of images and movie for all "2D"
    # files which contain grid stacks.  e.g. *.nc'
    #----------------------------------------------
    os.chdir( output_dir )
    nc_file_list = glob.glob('*.nc')
    for nc_file in nc_file_list:
        #------------------------------------------
        # Change the stretch for specific files ?
        #------------------------------------------
#         if nc_file.endswith('???.nc'):
#             cur_stretch = 'power3'
#             stretch = 'hist_equal'
          
        if ('IDs' in nc_file):
            LAND_SEA_BACKDROP     = True
            LAND_SEA_RED_BACKDROP = False
        else:
            LAND_SEA_BACKDROP     = False
            LAND_SEA_RED_BACKDROP = True
                      
        #------------------------------------
        # First, create a set of PNG images
        #-------------------------------------------------
        # Note: Linear stretch may be necessary to avoid
        #       having the land and sea change color.
        #       Power3 stretch may also be okay.
        #-------------------------------------------------
        save_grid_stack_as_images( nc_file, temp_png_dir,
                                   ### extent=None,  # auto-computed
                                   ## stretch='hist_equal', a=1, b=2, p=0.5,
                                   stretch='linear',
                                   ## stretch='power3', a=1, b=2, p=0.5,
                                   cmap='jet', REPORT=True,
                                   #### cmap='rainbow', REPORT=True, # jet > rainbow
                                   LAND_SEA_BACKDROP=LAND_SEA_BACKDROP,
                                   LAND_SEA_RED_BACKDROP=LAND_SEA_RED_BACKDROP,
                                   ## LAND_SEA_BACKDROP=True,  ####
                                   xsize=xsize2D, ysize=ysize2D, dpi=dpi)

        #----------------------------------------------
        # Create movie from set of images in temp_png
        #----------------------------------------------
        # movie_fps = "frames per second"
        mp4_file = nc_file.replace('.nc', '.mp4')
        mp4_path = (movie_dir + mp4_file)
        create_movie_from_images( mp4_path, temp_png_dir,
                                  fps=movie_fps, REPORT=True)

        #-----------------------------------
        # Delete all PNG files in temp_png
        #-----------------------------------
        delete_png_files( temp_png_dir )

    #---------------------------------------------------
    # Convert MP4 movies to WEBM movies with opacity ?
    # To hopefully allow backdrop overlay in Causemos.
    #---------------------------------------------------
    ### if (opacity == 1.0) or not(WEBM):
    if not(WEBM):
        os.chdir( start_dir )  # need this
        return

    #-----------------------------------------        
    # Get the full path to ffmpeg executable
    #--------------------------------------------
    # -y means OK to overwrite files w/o asking
    # -hide_banner suppresses extra chatter
    #--------------------------------------------   
    print('Converting mp4 to webm format with opacity...')
    ffmpeg_path = imff.get_ffmpeg_exe()
    hb = ' -hide_banner -y '
    
    #------------------------------------------------        
    # Convert mp4 to webm with opacity, with ffmpeg
    #------------------------------------------------
    os.chdir( movie_dir )
    mp4_file_list = glob.glob('*.mp4')
    for mp4_file in mp4_file_list:
        cmd_p1 = ffmpeg_path + hb + ' -i '
        cmd_p2 = ' -filter_complex "format=yuva444p,colorchannelmixer=aa=0.5" '
        webm_file = mp4_file.replace('.mp4', '.webm')
        # print('#### mp4_file  =', mp4_file)
        # print('#### webm_file =', webm_file)
        cmd       = (cmd_p1 + mp4_file + cmd_p2 + webm_file)
        os.system( cmd )
    print('Finished converting mp4 to webm format with opacity.')
    print('   opacity =', opacity )
    print()
    os.chdir( start_dir )  ####
    
#   create_media_files()
#----------------------------------------------------------------------------
def delete_png_files( temp_png_dir ):

    png_files = glob.glob( temp_png_dir + '*.png' )
    for file in png_files:
        try:
           os.remove( file )
        except OSError as e:
            print("Error: %s : %s" % (file, e.strerror)) 

    #-------------------
    # Is this needed ?
    #-------------------
    time.sleep( 1.0 )

    print('Finished deleting PNG files in:')
    print('  ' + temp_png_dir)
    print()

#   delete_png_files()
#----------------------------------------------------------------------------






   