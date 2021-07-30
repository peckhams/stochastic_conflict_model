#! /usr/bin/env python

import argparse
from conflict.model import conflict

def main(args=None):
    """The main routine."""

    parser = argparse.ArgumentParser(description='Run localized conflict model.')

    parser.add_argument('--cfg_file',
                        help='filename with path for the configuration (CFG) file')

    args = parser.parse_args()

    c = conflict.conflict()
    c.run_model( cfg_file=args.cfg_file )

#------------------------------------------
if __name__ == "__main__":
    main()