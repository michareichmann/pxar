#!/usr/bin/env python
# --------------------------------------------------------
#       Class to write a list is pXar Events into a root tree
# created on April 18th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TTree, TFile, vector
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from collections import OrderedDict
from os.path import isfile


class TreeWriter:

    def __init__(self, data):

        self.Data = data
        self.File = None
        self.Tree = None
        self.VectorBranches = self.init_vector_branches()

        self.RunFileName = 'runNumber.txt'
        self.RunNumber = self.load_run_number()

        self.ProgressBar = None

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()], maxval=n)
        self.ProgressBar.start()

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

    @staticmethod
    def init_vector_branches():
        dic = OrderedDict([('col', vector('unsigned short')()),
                           ('row', vector('unsigned short')()),
                           ('adc', vector('short')())])
        return dic

    def clear_vectors(self):
        for key in self.VectorBranches.iterkeys():
            self.VectorBranches[key].clear()

    def set_branches(self):
        for key, vec in self.VectorBranches.iteritems():
            self.Tree.Branch(key, vec)

    def write_tree(self):
        self.File = TFile('run{n}.root'.format(n=str(self.RunNumber).zfill(3)), 'RECREATE')
        self.Tree = TTree('tree', 'The source tree')
        self.set_branches()
        self.start_pbar(len(self.Data))
        for i, ev in enumerate(self.Data):
            self.ProgressBar.update(i + 1)
            self.clear_vectors()
            for pix in ev.pixels:
                self.VectorBranches['col'].push_back(int(pix.column))
                self.VectorBranches['row'].push_back(int(pix.row))
                self.VectorBranches['adc'].push_back(int(pix.value))
            self.Tree.Fill()
        self.ProgressBar.finish()
        self.File.cd()
        self.File.Write()
        self.File.Close()
        self.save_run_number()


if __name__ == '__main__':
    z = TreeWriter(None)
