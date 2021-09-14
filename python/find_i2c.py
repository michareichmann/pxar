#!/usr/bin/env python3
# --------------------------------------------------------
#       Script to find the I2C of a ROC
# created on February 20th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from argparse import ArgumentParser
from helpers.files import find_i2c


p = ArgumentParser()
p.add_argument('n', nargs='?', default=16, type=int)
args = p.parse_args()

find_i2c(args.n)
