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
        branches = []
        for i in xrange(self.NPlanes):
            branches.append(OrderedDict([('Value', zeros(50, 'f8')),
                                         ('Timing', zeros(50, 'f8')),
                                         ('PixX', zeros(50, 'i4')),
                                         ('PixY', zeros(50, 'i4')),
                                         ('HitInCluster', zeros(50, 'i4')),
                                         ('TriggerCount', zeros(50, 'i4'))]))
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
        self.EventBranches['TimeStamp'][0] = time() * 1000
        n_hits = zeros(self.NPlanes)
        x, y, adc = [[]] * self.NPlanes, [[]] * self.NPlanes, [[]] * self.NPlanes
        for pix in ev.pixels:
            n_hits[pix.roc] += 1
            x[pix.roc].append(pix.column)
            y[pix.roc].append(pix.row)
            adc[pix.roc].append(pix.value)
        for i in xrange(self.NPlanes):
            self.VectorBranches[i]['PixX'] = array(x[i], 'i4')
            self.VectorBranches[i]['PixY'] = array(y[i], 'i4')
            self.VectorBranches[i]['Value'] = array(adc[i], 'f8')
        for j in xrange(self.NPlanes):
            self.ScalarBranches[j]['NHits'][0] = n_hits[j]
            self.VectorBranches[j]['Timing'] = array([ev.triggerPhases[i] for i in xrange(len(ev.header))], 'f8')
            self.VectorBranches[j]['TriggerCount'] = array([ev.triggerCounts[i] for i in xrange(len(ev.header))], 'i4')
        for i in xrange(self.NPlanes):
            self.Trees[i].Fill()
        self.EventTree.Fill()
        self.NEvents += 1


if __name__ == '__main__':
    z = TreeWriter(None)
