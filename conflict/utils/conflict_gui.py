"""
This module defines a class called "conflict_gui" that can be used to
create a graphical user interface (GUI) for the Localized Conflict
Model within a Jupyter notebook.  If used with Binder, this GUI runs
in a browser window and does not require the user to install anything
on their computer.  However, this module should be included in the
same directory as the Jupyter notebook.
"""
#------------------------------------------------------------------------
#
#  Copyright (C) 2021.  Scott D. Peckham
#
#------------------------------------------------------------------------

# Moved conflict.py into conflict/model
# since package and module have same name.
from conflict.model import conflict

import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display, HTML
## from IPython.core.display import display
## from IPython.lib.display import display

import numpy as np
import os, os.path
# import conflict_plot as cp

#------------------------------------------------------------------------
#
#  class conflict_gui
#      __init__()
#      pix_str()
#      show_gui()
#      make_tab_gui()
#      get_padding()
#      make_input_panel()
#      reset_input_panel()
#      make_output_panel()
#
#------------------------------------------------------------------------
class conflict_gui:
    #--------------------------------------------------------------------
    def __init__(self):

        self.version  = '0.5'
        self.user_var = None
        #----------------------------------------------------------
        # "full_box_width" = (label_width + widget_width)
        # gui_width = left_label_width + mid_width + button_width 
        # The 2nd, label + widget box, is referred to as "next".
        # (2 * half_widget_width) + left_label + next_label = 540
        #----------------------------------------------------------       
        # self.gui_width         = 680  # (orig)
        self.gui_width         = 600
        self.left_label_width  = 150
        self.next_label_width  = 70
        self.all_label_width   = 170
        # self.full_box_width    = 540  # (orig)
        self.full_box_width    = 550
        self.widget_width      = (self.full_box_width - self.left_label_width)
        # self.half_widget_width = (self.full_box_width - self.all_label_width)/2
        # self.half_widget_width = 183
        self.left_widget_width = 150
        self.next_widget_width = 150
        self.left_box_width    = (self.left_label_width + self.left_widget_width)
        self.next_box_width    = (self.next_label_width + self.next_widget_width)
        self.button_width      = 70   # big enough for "Reset"
        #-----------------------------------------------------
        self.default_U_file    = 'input_files/Horn_of_Africa_GPW-v4_pop_count_2020_450sec.tif'
        self.default_C1_file   = '(none, uniform)'
        self.default_C2_file   = '(none, uniform)'
        self.default_gui_cfg_file = '~/Conflict/Input/gui_conflict.cfg'
        self.default_out_file  = '~/Conflict/Output/conflicts.rts'
        self.default_IDs_file  = '~/Conflict/Output/conflict_IDs.rts'
        #-----------------------------------------------------
        self.gui_width_px  = self.pix_str( self.gui_width )
        #---------------------------------------------------
        # These styles are used to control width of labels
        # self.init_label_style is the initial default.
        #---------------------------------------------------
        llw_px = self.pix_str( self.left_label_width )
        nlw_px = self.pix_str( self.next_label_width )
        self.init_label_style = {'description_width': 'initial'}
        self.left_label_style = {'description_width': llw_px}
        self.next_label_style = {'description_width': nlw_px}
        
    #   __init__()
    #--------------------------------------------------------------------
    def pix_str(self, num):
        return str(num) + 'px'

    #--------------------------------------------------------------------
    def show_gui(self):

        #------------------------------------   
        # Create & display the complete GUI
        #-----------------------------------
        self.make_tab_gui()

        gui_output = widgets.Output()
        display(self.gui, gui_output)
  
    #   show_gui()
    #--------------------------------------------------------------------
    def make_tab_gui(self):

        gui_width_px = self.gui_width_px
        
        self.make_input_panel()
        self.make_output_panel()
        
#         self.make_map_panel( SHOW_MAP=SHOW_MAP )
#         self.make_datetime_panel()
#         self.make_download_panel()
#         self.make_prefs_panel()
        #---------------------------
        p0 = self.input_panel
        p1 = self.output_panel
#         p1 = self.map_panel
#         p2 = self.datetime_panel
#         p3 = self.download_panel
#         p4 = self.prefs_panel
        #---------------------------
        p0_title = 'Inputs'
        p1_title = 'Outputs'
#         p0_title = 'Browse Data'
#         p1_title = 'Spatial Extent'
#         p2_title = 'Date Range'
#         p3_title = 'Download Data'
#         p4_title = 'Settings'
          
        #-------------------------------------------------------
        # selected_index=0 shows Browse Data panel
        #-------------------------------------------------------
        # tab = widgets.Tab( children=[p0, p1, p2, p3, p4], 
        tab = widgets.Tab( children=[p0, p1],        
                           selected_index=0,
                           layout=Layout(width=gui_width_px) )
        tab.set_title(0, p0_title)
        tab.set_title(1, p1_title)
#         tab.set_title(2, p2_title)
#         tab.set_title(3, p3_title)
#         tab.set_title(4, p4_title)
        #### tab.titles = [str(i) for i in range(len(children))]
        
        # L_tags = "<b><font size=5>"
        # R_tags = "</font></b>"
        # heading = (L_tags + title + R_tags)
        pad  = self.get_padding(1, HORIZONTAL=False)  # 1 lines
        head = widgets.HTML(value=f"<font size=5>Stochastic Conflict Model User Interface</font>")
        # head = widgets.Label('Conflict 0.5 User Interface')
        ## self.gui = widgets.VBox([pad, head, acc])
        self.gui = widgets.VBox([head, tab])   # (no padding above)

    #   make_tab_gui()
    #--------------------------------------------------------------------
    def get_padding(self, n, HORIZONTAL=True):
 
        #-------------------------------       
        # Get some white space padding
        #-------------------------------
        if (HORIZONTAL):
            #--------------------------------
            # Use overloaded multiplication
            #--------------------------------
            ## s  = (' ' * n)  # overloaded multiplication
            s  = "<p>" + ('&nbsp;' * n) + "</p>"
            pad = widgets.HTML( value=s )
        else:
            s = ("<br>" * n)
            pad = widgets.HTML( value=s )
        return pad
        
    #   get_padding()
    #--------------------------------------------------------------------
    def make_input_panel(self):

        #-----------------------------------
        # Browse data on an OpenDAP server
        #-----------------------------------
        left_style     = self.left_label_style
        next_style     = self.next_label_style
        full_width_px  = self.pix_str( self.full_box_width )
        left_width_px  = self.pix_str( self.left_box_width )
        next_width_px  = self.pix_str( self.next_box_width )
        btn_width_px   = self.pix_str( self.button_width )
        #---------------------------------------------------------------
        o1 = widgets.Text(description='N_time_steps:', style=left_style,
                          value='100', layout=Layout(width=left_width_px) )                  
#         o2 = widgets.Text(description='Grid ncols:', style=next_style,
#                           value='', layout=Layout(width=next_width_px) )
        o2 = widgets.Text(description='Grid ncols:', style=left_style,
                          value='240', layout=Layout(width=left_width_px) )        
        o3 = widgets.Text(description='Grid nrows:', style=next_style,
                          value='240', layout=Layout(width=next_width_px) ) 
        #---------------------------------------------------------------
        o4 = widgets.Text(description='Unrest Grid File:',
                          value=self.default_U_file,  #########
                          disabled=False, style=left_style,
                          layout=Layout(width=full_width_px))
        # b1 = widgets.Button(description="Go", layout=Layout(width=btn_width_px))
        o5 = widgets.Text(description='Connectivity File 1:',
                          value=self.default_C1_file,  #########
                          disabled=False, style=left_style,
                          layout=Layout(width=full_width_px))
        o6 = widgets.Text(description='Connectivity File 2:',
                          value=self.default_C2_file,  #########
                          disabled=False, style=left_style,
                          layout=Layout(width=full_width_px))
        #---------------------------------------------------------------
        o7 = widgets.FloatSlider(description='Emergence factor:', value=0.005,
                                 min=0.001, max=1.0, step=0.001,
                                 disabled=False,
                                 continuous_update=False,
                                 orientation='horizontal',
                                 readout=True, readout_format='.3f',
                                 style=left_style,
                                 layout=Layout(width=full_width_px))
        o8 = widgets.FloatSlider(description='Spreading factor:', value=0.2,
                                 min=0.0, max=1.0, step=0.001,
                                 disabled=False,
                                 continuous_update=False,
                                 orientation='horizontal',
                                 readout=True, readout_format='.3f',
                                 style=left_style,
                                 layout=Layout(width=full_width_px))
        o9 = widgets.FloatSlider(description='Resolution probability:', value=0.47,
                                 min=0.0, max=1.0, step=0.001,
                                 disabled=False,
                                 continuous_update=False,
                                 orientation='horizontal',
                                 readout=True, readout_format='.3f',
                                 style=left_style,
                                 layout=Layout(width=full_width_px))                                              
        o10 = widgets.Text(description='Status:', style=left_style,
                          value='Ready.', layout=Layout(width=full_width_px) )
        b2 = widgets.Button(description="Run", layout=Layout(width=btn_width_px))            
        ## b2 = widgets.Button(description="Reset", layout=Layout(width=btn_width_px))
        ## pd = widgets.HTML(('&nbsp;' * 1))  # for padding
        
        #-------------------------------
        # Arrange widgets in the panel
        #-------------------------------
        dims_box = widgets.HBox([o2, o3])
        ## U_file_box = widgets.HBox([o1, b1])      # U_file + Go button
        stat_box = widgets.HBox([o10, b2])      # status + Run button
        ## pad_box  = widgets.VBox([pd, pd])
        ## mid_box  = widgets.HBox([name_box, unit_box])
        ## mid_box  = widgets.HBox([name_box, pad_box, unit_box])
        panel = widgets.VBox([o1, dims_box, o4, o5, o6, o7, o8, o9, stat_box])
 
        self.input_n_steps   = o1
        self.input_n_cols    = o2
        self.input_n_rows    = o3
        self.input_U_file    = o4
        self.input_C1_file   = o5
        self.input_C2_file   = o6
        self.input_c_emerge  = o7
        self.input_c_spread  = o8
        self.input_p_resolve = o9
        self.input_status    = o10
        self.input_panel     = panel

        #-----------------
        # Event handlers
        #-----------------------------------------------------
        # Note: NEED to set names='value' here.  If names
        #       keyword is omitted, only works intermittently.
        #------------------------------------------------------------
        # "on_click" handler function is passed b1 as argument.
        # "observe" handler function is passed "change", which
        # is a dictionary, as argument. See Traitlet events.
        #------------------------------------------------------------
        b2.on_click( self.run_conflict_model )  ###########
        # b1.on_click( self.reset_input_panel )
        # o2.observe( self.update_data_panel, names=['options','value'] )
        # o3.observe( self.update_var_info, names=['options', 'value'] )
        ## o3.observe( self.update_var_info, names='value' )
        ## o2.observe( self.update_data_panel, names='All' )
        ## o3.observe( self.update_var_info, names='All' )
 
        #-------------------------------------------------------    
        # It turned out this wasn't an issue, but interesting.
        #-------------------------------------------------------
        # Note: Method functions have type "method" instead
        #       of "function" and therefore can't be passed
        #       directly to widget handlers like "on_click".
        #       But we can use the "__func__" attribute.
        #-------------------------------------------------------           
#         b1.on_click( self.update_filename_list.__func__ )
#         o2.observe( self.update_data_panel.__func__ )
#         o3.observe( self.update_var_info.__func__, names='value' )
       
    #   make_input_panel()
    #--------------------------------------------------------------------
    def reset_input_panel(self, caller_obj=None):

        #----------------------------------------------------
        # Note: This is called by the "on_click" method of
        # the "Reset" button beside the status box.
        # In this case, type(caller_obj) =
        # <class 'ipywidgets.widgets.widget_button.Button'>
        #----------------------------------------------------
#         self.input_U_file 
#         self.data_filename.options    = ['']
#         self.data_var_name.options    = ['']  # short names
#         self.data_var_long_name.value = ''
#         self.data_var_units.value     = ''
#         self.data_var_shape.value     = ''
#         self.data_var_dims.value      = ''
#         self.data_var_type.value      = ''
#         self.data_var_atts.options    = ['']
        self.input_status.value        = 'Ready.'   
        #------------------------------------------
        # self.download_log.value = ''

    #   reset_input_panel() 
    #--------------------------------------------------------------------
    def make_output_panel(self):

        #-----------------------------------
        # Browse data on an OpenDAP server
        #-----------------------------------
        left_style     = self.left_label_style
        next_style     = self.next_label_style
        full_width_px  = self.pix_str( self.full_box_width )
        left_width_px  = self.pix_str( self.left_box_width )
        next_width_px  = self.pix_str( self.next_box_width )
        btn_width_px   = self.pix_str( self.button_width )
        #-------------------------------------------------------------
        o1 = widgets.Text(description='GUI-created cfg file:',
                          value=self.default_gui_cfg_file,
                          disabled=False, style=left_style,
                          layout=Layout(width=full_width_px))
        o2 = widgets.Text(description='Conflict RTS File:',
                          value=self.default_out_file,
                          disabled=False, style=left_style,
                          layout=Layout(width=full_width_px))
        o3 = widgets.Text(description='Conflict IDs RTS File:',
                          value=self.default_IDs_file,
                          disabled=False, style=left_style,
                          layout=Layout(width=full_width_px))
        
        #-------------------------------
        # Arrange widgets in the panel
        #-------------------------------
        panel = widgets.VBox([o1, o2, o3])
 
        self.output_gui_cfg_file  = o1
        self.output_conflict_file = o2
        self.output_IDs_file      = o3 
        self.output_panel         = panel

        #-----------------
        # Event handlers
        #-----------------------------------------------------
        # Note: NEED to set names='value' here.  If names
        #       keyword is omitted, only works intermittently.
        #------------------------------------------------------------
        # "on_click" handler function is passed b1 as argument.
        # "observe" handler function is passed "change", which
        # is a dictionary, as argument. See Traitlet events.
        #------------------------------------------------------------
        # b2.on_click( self.run_conflict_model )  ###########
        # b1.on_click( self.reset_input_panel )
        # o2.observe( self.update_data_panel, names=['options','value'] )
        # o3.observe( self.update_var_info, names=['options', 'value'] )
       
    #   make_output_panel()
    #-------------------------------------------------------------------- 
    def write_config_file(self):
    
        #-----------------------------
        # Open new cfg file to write
        #-----------------------------
        cfg_file = self.output_gui_cfg_file.value
        cfg_file = os.path.expanduser( cfg_file )
        
        cfg_unit = open( cfg_file, 'w')
        cfg_unit.write('#-----------------------------------------------\n')
        cfg_unit.write('# Configuration file for Conflict Model v. 0.5\n')
        cfg_unit.write('#-----------------------------------------------\n')
        cfg_unit.write('n_steps   = ' + str(self.input_n_steps.value) + '\n')
        cfg_unit.write('n_cols    = ' + str(self.input_n_cols.value) + '\n')
        cfg_unit.write('n_rows    = ' + str(self.input_n_rows.value) + '\n')
        #--------------------------------------------------------------------
        U_file = str(self.input_U_file.value)
        if (U_file == '(none, uniform)'):
            cfg_unit.write("U_file    = ''\n")
        else:
            cfg_unit.write('U_file    = ' + "'" + U_file + "'\n")
        #--------------------------------------------------------------------
        C1_file = str(self.input_C1_file.value)
        if (C1_file == '(none, uniform)'):
            cfg_unit.write("C1_file   = ''\n")
        else:           
            cfg_unit.write('C1_file   = ' + "'" + C1_file + "'\n")
        #--------------------------------------------------------------------        
        C2_file = str(self.input_C2_file.value)
        if (C2_file == '(none, uniform)'):
            cfg_unit.write("C2_file   = ''\n")
        else:           
            cfg_unit.write('C2_file   = ' + "'" + C2_file + "'\n") 
        #--------------------------------------------------------------------        
        cfg_unit.write('c_emerge  = ' + str(self.input_c_emerge.value) + '\n')
        cfg_unit.write('c_spread  = ' + str(self.input_c_spread.value) + '\n')     
        cfg_unit.write('p_resolve = ' + str(self.input_p_resolve.value) + '\n')
        cfg_unit.write('spread_method = ' + str(1) + '\n')
        cfg_unit.write('time_lag  = ' + str(1) + '\n')
        cfg_unit.write('REPORT    = ' + str(1) + '\n')                       
        cfg_unit.write("out_file  = '" + self.output_conflict_file.value + "'\n")
        cfg_unit.write("IDs_file  = '" + self.output_IDs_file.value + "'\n")         
        cfg_unit.close()
         
        self.gui_cfg_file = cfg_file

    #   write_config_file()
    #-------------------------------------------------------------------- 
    def run_conflict_model(self, caller_obj=None ): 

        #----------------------------------------------
        # Note:  Need the "caller_obj" argument here.
        #----------------------------------------------
        self.input_status.value = 'Writing new cfg file...'  
        self.write_config_file()
        
        cfg_file = self.gui_cfg_file
        self.input_status.value = 'Instantiating conflict model...'
        c = conflict.conflict()
        self.input_status.value = 'Running conflict model...'         
        c.run_model( cfg_file=cfg_file )
        self.input_status.value = 'Finished.' 
        
    #   run_conflict_model()
    #--------------------------------------------------------------------
    
        