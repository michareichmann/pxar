#!/usr/bin/env python
# --------------------------------------------------------
#       Class to read pXar events from file
# created on November 13th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from glob import glob
from helpers.draw import *


class FileReader:

    def __init__(self, run_number, config_name, file_type):

        self.Dir = dirname(dirname(realpath(__file__)))
        self.Config = load_config(join(self.Dir, 'config', config_name))
        self.DataDir = self.Config.get('MAIN', 'data directory')
        self.RunNumber = self.load_run_number(run_number)
        self.FileName = '{}_{:03d}.{}'.format(self.Config.get('MAIN', 'filename'), self.RunNumber, file_type)
        self.NPlanes = self.Config.getint('MAIN', 'number of planes')

        self.NCols = self.Config.getint('CHIP', 'columns')
        self.NRows = self.Config.getint('CHIP', 'rows')

        self.PBar = PBar()

        self.File = self.load_file()
        self.Plotter = Draw()

    def load_run_number(self, run_number):
        return int(remove_letters(max(glob(join(self.DataDir, '*'))).split('.')[0])) if run_number is None else int(run_number)

    def load_file(self):
        pass
