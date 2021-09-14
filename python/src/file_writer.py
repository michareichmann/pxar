#!/usr/bin/env python
# --------------------------------------------------------
#       Class to save pXar events to file
# created on November 12th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from helpers.utils import *
from os.path import join, dirname, realpath
from glob import glob
from ROOT import TF1
from numpy import arange, split, genfromtxt, array


class FileWriter:

    def __init__(self, config_name, file_type):

        self.Dir = dirname(dirname(realpath(__file__)))
        self.Config = load_config(join(self.Dir, 'config', config_name))
        self.DataDir = self.Config.get('MAIN', 'data directory')
        self.RunNumber = self.load_run_number()
        self.FileName = '{}_{:03d}.{}'.format(self.Config.get('MAIN', 'filename'), self.RunNumber, file_type)
        self.NPlanes = self.Config.getint('MAIN', 'number of planes')

        self.WBC = self.Config.getint('CHIP', 'wbc')
        self.NCols = self.Config.getint('CHIP', 'columns')
        self.NRows = self.Config.getint('CHIP', 'rows')
        self.Trim = self.Config.get('CHIP', 'trim') if self.Config.has_option('CHIP', 'trim') else ''
        self.NEvents = 0

        # Pulse Height Calibrations
        self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
        self.Parameters = self.load_calibration_fitpars()

        self.PBar = PBar()

    def __del__(self):
        self.save_file()

    def load_run_number(self):
        f_names = glob(join(self.DataDir, '*'))
        return 0 if not f_names else max(int(remove_letters(f_name.split('.')[0])) for f_name in f_names) + 1

    def load_calibration_fitpars(self):
        split_at = arange(self.NRows, self.NCols * self.NRows, self.NRows)  # split at every new column (after n_rows)
        lst = array([split(genfromtxt(filename, skip_header=3, usecols=arange(4)), split_at) for filename in glob('phCalibrationFitErr{}*'.format(self.Trim))])
        if not lst.size or not lst[0].size:
            warning('Did not find calibration file for trim {}! '.format(self.Trim))
            return
        return lst

    def get_vcal(self, roc, col, row, adc):
        if self.Parameters is None:
            return adc
        self.Fit.SetParameters(*self.Parameters[roc][col][row])
        return self.Fit.GetX(adc)

    def save_file(self):
        pass
