#!/usr/bin/env python
# --------------------------------------------------------
#       Class to write a list is pXar Events into a root tree
# created on February 20th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import vector
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

        self.HitDirs = []
        self.Trees = []
        self.init_trees()

        self.set_event_branches()
        self.set_branches()

    @staticmethod
    def init_event_branches():
        return OrderedDict([('TimeStamp', array([0], 'i8'))])

    def init_scalar_branches(self):
        return [OrderedDict([('NHits', array([0], 'i'))]) for _ in xrange(self.NPlanes)]

    def init_vector_branches(self):
        branches = []
        for i in xrange(self.NPlanes):
            branches.append(OrderedDict([('Value', vector('double')()),
                                         ('Timing', vector('double')()),
                                         ('PixX', vector('int')()),
                                         ('PixY', vector('int')()),
                                         ('HitInCluster', vector('int')()),
                                         ('TriggerCount', vector('int')())]))
        return branches

    def set_event_branches(self):
        for key, value in self.EventBranches.iteritems():
            print '{}/{}'.format(key, type_dict[value[0].dtype.name])
            self.EventTree.Branch(key, value, '{}/{}'.format(key, type_dict[value[0].dtype.name]))

    def init_trees(self):
        for i in xrange(self.NPlanes):
            self.HitDirs.append(self.File.mkdir('Plane{}'.format(i + self.Config.getint('TREE', 'plane number'))))
            self.HitDirs[i].cd()
            self.Trees.append(self.init_tree())


    def write(self, ev):
        self.clear_vectors()
        n_hits = int(len(ev.pixels))
        self.ScalarBranches[ev.roc]['NHits'][0] = n_hits
        self.EventBranches['TimeStamp'][0] = long(time() * 1000)
        for pix in ev.pixels:
            self.VectorBranches[ev.roc]['PixX'].push_back(int(pix.column))
            self.VectorBranches[ev.roc]['PixY'].push_back(int(pix.row))
            self.VectorBranches[ev.roc]['Value'].push_back(int(pix.value))
        for i in xrange(len(ev.header)):
            self.VectorBranches[ev.roc]['Timing'].push_back(ev.triggerPhases[i])
            self.VectorBranches[ev.roc]['TriggerCount'].push_back(ev.triggerCounts[i])
        for i in xrange(self.NPlanes):
            self.Trees[i].Fill()
        self.EventTree.Fill()


if __name__ == '__main__':
    z = TreeWriter(None)
