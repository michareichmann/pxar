#!/usr/bin/env python3
# --------------------------------------------------------
#       Script to change the I2C of a ROC
# created on February 20th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.files import rename_files, change_config, get_old_i2c
from argparse import ArgumentParser


p = ArgumentParser(usage='change_i2c new_i2c, change_i2c old_i2c new_i2c')
p.add_argument('n')
p.add_argument('o', nargs='?', default=None, help='default: smallest i2c')
args = p.parse_args()

n, o = (args.n, get_old_i2c()) if args.o is None else (args.o, get_old_i2c(args.n))
rename_files(n, o)
change_config(n, o)
