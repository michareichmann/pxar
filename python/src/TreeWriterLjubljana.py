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
        return OrderedDict([('TimeStamp', array([0], 'f8'))])

    def init_scalar_branches(self):
        return [OrderedDict([('NHits', array([0], 'i'))]) for _ in xrange(self.NPlanes)]

    def init_vector_branches(self):
        n = 5000
        branches = []
        for i in xrange(self.NPlanes):
            branches.append(OrderedDict([('Value', zeros(n, 'f8')),
                                         ('Timing', zeros(n, 'f8')),
                                         ('PixX', zeros(n, 'i4')),
                                         ('PixY', zeros(n, 'i4')),
                                         ('HitInCluster', zeros(n, 'i4')),
                                         ('TriggerCount', zeros(n, 'i4'))]))
        return branches

    def set_event_branches(self):
        for key, value in self.EventBranches.iteritems():
            self.EventTree.Branch(key, value, '{}/{}'.format(key, type_dict[value[0].dtype.name]))

    def init_trees(self):
        for i in xrange(self.NPlanes):
            self.File.cd()
            self.HitDirs.append(self.File.mkdir('Plane{}'.format(i + self.Config.getint('TREE', 'plane number'))))
            self.HitDirs[i].cd()
            self.Trees.append(self.init_tree())

    def set_branches(self):
        for itree in xrange(self.NPlanes):
            for key, value in self.ScalarBranches[itree].iteritems():
                self.Trees[itree].Branch(key, value, '{}/{}'.format(key, type_dict[value[0].dtype.name]))
            for key, vec in self.VectorBranches[itree].iteritems():
                array_size = '1' if key in ['Timing', 'TriggerCount'] else 'NHits'
                self.Trees[itree].Branch(key, vec, '{key}[{n}]/{type}'.format(key=key, n=array_size, type=type_dict[vec[0].dtype.name]))

    def write(self, ev):
        self.EventBranches['TimeStamp'][0] = time() * 1000
        n_hits = zeros(self.NPlanes, 'i')
        x, y, adc = [[[] for _ in xrange(self.NPlanes)] for _ in xrange(3)]
        for pix in ev.pixels:
            n_hits[pix.roc] += 1
            x[pix.roc].append(pix.column)
            y[pix.roc].append(pix.row)
            adc[pix.roc].append(pix.value)
        for i in xrange(self.NPlanes):
            self.ScalarBranches[i]['NHits'][0] = n_hits[i]
            for j in xrange(len(ev.header)):
                self.VectorBranches[i]['Timing'][j] = ev.triggerPhases[j]
                self.VectorBranches[i]['TriggerCount'][j] = ev.triggerCounts[j]
            for j in xrange(n_hits[i]):
                self.VectorBranches[i]['PixX'][j] = x[i][j]
                self.VectorBranches[i]['PixY'][j] = y[i][j]
                self.VectorBranches[i]['Value'][j] = adc[i][j]
            self.Trees[i].Fill()
        self.EventTree.Fill()
        self.NEvents += 1


if __name__ == '__main__':
    z = TreeWriter(None)
