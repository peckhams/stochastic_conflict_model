#
# Copyright (c) 2021, Scott D. Peckham
#
# run_model()
#
from conflict.model import conflict

#-----------------------------------------------------------------------
def run_model( cfg_file, SILENT=False):

    #----------------------------------------------------------
    # Note: 
    #----------------------------------------------------------
    c = conflict.conflict()
    c.run_model( cfg_file=cfg_file )

#   run_model()
#-------------------------------------------------------------------

