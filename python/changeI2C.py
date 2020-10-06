#!/usr/bin/env python
from glob import glob
from shutil import move
from argparse import ArgumentParser

p = ArgumentParser()
p.add_argument('n')

args = p.parse_args()

o = ''
for name in glob('*_C*'):
    o = name[name.index('_C') + 2]
    move(name, name.replace('_C' + o, '_C' + args.n))

f = open('configParameters.dat', 'r+')
lines = []
for line in f.readlines():
    if line.startswith('nRocs'):
        lines.append(line.replace('i2c: ' + o, 'i2c: ' + args.n))
    else:
        lines.append(line)

f.seek(0)
f.writelines(lines)
f.truncate()
f.close()
