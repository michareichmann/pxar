#!/usr/bin/env python
# --------------------------------------------------------
#       Class to write a list is pXar Events into a root tree
# created on February 20th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import vector, gRandom
from collections import OrderedDict
from TreeWriter import *
from time import time
from numpy import array


class TreeWriterLjubljana(TreeWriter):

    def __init__(self, config_name='ljutel'):

        TreeWriter.__init__(self, config_name)

        self.File.cd()
        self.EventTree = TTree('Event', 'Event Information')
        self.EventBranches = self.init_event_branches()

        self.HitDir = self.File.mkdir('Plane{}'.format(self.Config.get('TREE', 'plane number')))
        self.HitDir.cd()
        self.Tree = self.init_tree()

        self.set_event_branches()
        self.set_branches()

    @staticmethod
    def init_event_branches():
        return OrderedDict([('TimeStamp', array([0], 'i8'))])

    @staticmethod
    def init_scalar_branches():
        return OrderedDict([('NHits', array([0], 'i'))])

    @staticmethod
    def init_vector_branches():
        return OrderedDict([('Value', vector('double')()),
                            ('Timing', vector('double')()),
                            ('PixX', vector('int')()),
                            ('PixY', vector('int')()),
                            ('HitInCluster', vector('int')()),
                            ('TriggerCount', vector('int')())])

    def set_event_branches(self):
        for key, value in self.EventBranches.iteritems():
            print '{}/{}'.format(key, type_dict[value[0].dtype.name])
            self.EventTree.Branch(key, value, '{}/{}'.format(key, type_dict[value[0].dtype.name]))


    def write(self, ev):
        self.clear_vectors()
        n_hits = int(len(ev.pixels))
        self.ScalarBranches['NHits'][0] = n_hits
        self.EventBranches['TimeStamp'][0] = long(time() * 1000)
        for pix in ev.pixels:
            self.VectorBranches['PixX'].push_back(int(pix.column))
            self.VectorBranches['PixY'].push_back(int(pix.row))
            self.VectorBranches['Value'].push_back(int(pix.value))
        for i in xrange(len(ev.header)):
            self.VectorBranches['Timing'].push_back(ev.triggerPhases[i])
            self.VectorBranches['TriggerCount'].push_back(ev.triggerCounts[i])
        self.Tree.Fill()
        self.EventTree.Fill()


if __name__ == '__main__':
    z = TreeWriter(None)
