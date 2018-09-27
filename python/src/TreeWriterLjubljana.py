#!/usr/bin/env python
# --------------------------------------------------------
#       Class to write a list is pXar Events into a root tree
# created on February 20th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import vector
from collections import OrderedDict
from TreeWriter import *
from time import time
from numpy import array, zeros


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
            self.File.cd()
            self.HitDirs.append(self.File.mkdir('Plane{}'.format(i + self.Config.getint('TREE', 'plane number'))))
            self.HitDirs[i].cd()
            self.Trees.append(self.init_tree())


    def write(self, ev):
        self.clear_vectors()
        self.EventBranches['TimeStamp'][0] = long(time() * 1000)
        n_hits = zeros(self.NPlanes)
        for pix in ev.pixels:
            n_hits[pix.roc] += 1
            self.VectorBranches[pix.roc]['PixX'].push_back(int(pix.column))
            self.VectorBranches[pix.roc]['PixY'].push_back(int(pix.row))
            self.VectorBranches[pix.roc]['Value'].push_back(int(pix.value))
        for j in xrange(self.NPlanes):
            self.ScalarBranches[j]['NHits'][0] = n_hits[j]
            for i in xrange(len(ev.header)):
                self.VectorBranches[i]['Timing'].push_back(ev.triggerPhases[i])
                self.VectorBranches[i]['TriggerCount'].push_back(ev.triggerCounts[i])
        for i in xrange(self.NPlanes):
            self.Trees[i].Fill()
        self.EventTree.Fill()
	self.NEvents += 1


if __name__ == '__main__':
    z = TreeWriter(None)
