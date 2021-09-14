#!/usr/bin/env python
# --------------------------------------------------------
#       Class to write a list is pXar Events into a root tree
# created on February 20th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TTree, TFile
from helpers.utils import *
from os.path import join, dirname, realpath, basename
from time import sleep


class TreeWriter:

    def __init__(self, config_name):

        self.Dir = dirname(dirname(realpath(__file__)))
        self.Config = load_config(join(self.Dir, 'config', config_name))
        self.DataDir = self.Config.get('MAIN', 'data directory')
        self.RunFileName = join(self.DataDir, self.Config.get('MAIN', 'run number file'))
        self.RunNumber = self.load_run_number()
        self.FileName = self.Config.get('MAIN', 'filename')

        self.NPlanes = self.Config.getint('TREE', 'number of planes')

        # init Trees and Branches
        self.File = self.init_file()
        self.Trees = None  # the time of creation determines the directory structure...
        self.VectorBranches = self.init_vector_branches()
        self.ScalarBranches = self.init_scalar_branches()

        self.NEvents = 0

    def __del__(self):
        self.File.cd()
        self.File.Write()
        self.File.Close()
        self.save_run_number()
        info('Successfully saved the tree: {}'.format(basename(self.File.GetName())))
        sleep(.1)
        self.copy_file()

    def load_run_number(self):
        if isfile(self.RunFileName):
            f = open(self.RunFileName)
            run_number = int(f.readline())
            f.close()
            return run_number
        else:
            f = open(self.RunFileName, 'w')
            f.write('1')
            f.close()
            return 1

    def save_run_number(self):
        f = open(self.RunFileName, 'w')
        f.write('{n}'.format(n=self.RunNumber + 1))
        f.close()

    def init_file(self):
        ensure_dir(self.DataDir)
        return TFile(join(self.DataDir, '{}_{}.root'.format(self.FileName, str(self.RunNumber).zfill(3))), 'RECREATE')

    def init_tree(self):
        return TTree(*[self.Config.get('TREE', 'name')] * 2)

    @staticmethod
    def init_scalar_branches():
        return {}

    @staticmethod
    def init_vector_branches():
        return {}

    def clear_vectors(self):
        for i in range(self.NPlanes):
            for key in self.VectorBranches[i].iterkeys():
                self.VectorBranches[i][key].clear()

    def set_branches(self):
        for itree in range(self.NPlanes):
            for key, value in self.ScalarBranches[itree].iteritems():
                self.Trees[itree].Branch(key, value, '{}/{}'.format(key, type_dict[value[0].dtype.name]))
            for key, vec in self.VectorBranches[itree].iteritems():
                self.Trees[itree].Branch(key, vec, '{}[NHits]/{}'.format(key, type_dict[vec[0].dtype.name]))

    def write(self, event):
        pass

    def copy_file(self):
        pass


if __name__ == '__main__':
    z = TreeWriter(None)
