#!/usr/bin/env python3
# --------------------------------------------------------
#       small script to clean pXar data directory
# created on June 4th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.files import remove_logs, remove_swaps, remove_trim_files
from argparse import ArgumentParser


p = ArgumentParser()
p.add_argument('--all', '-a', action='store_true')
p.add_argument('-i2c', nargs='?', default=None)
args = p.parse_args()

remove_logs()
remove_swaps()
if args.all:
    remove_trim_files(args.i2c)
