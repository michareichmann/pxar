#!/usr/bin/env python
# --------------------------------------------------------
#       Class to save pXar events to hdf5 file
# created on November 12th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

import h5py
from numpy import cumsum, mean, sum, empty
from src.file_writer import *


class HDF5Writer(FileWriter):

    def __init__(self, config_name):
        FileWriter.__init__(self, config_name, 'hdf5')

        # Data
        self.NHits = []
        self.Hits = self.init_list()
        self.Clusters = self.init_list()
        self.NClusters = self.init_list()
        self.TriggerPhase = []

    # noinspection PyTypeChecker
    def init_list(self):
        v = empty(self.NPlanes, list)
        v.fill([])
        return v

    def convert(self):
        self.make_arrays()
        self.clusterise()

    def save_file(self):
        if len(self.NHits):
            ensure_dir(self.DataDir)
            info('saving file: {}'.format(self.FileName))
            with h5py.File(join(self.DataDir, self.FileName), 'w') as f:
                f.create_dataset('trigger_phase', data=array(self.TriggerPhase, 'u1')) if self.TriggerPhase else do_nothing()
                for roc in range(self.NPlanes):
                    grp = f.create_group('ROC{}'.format(roc))
                    grp.create_dataset('hits', data=self.Hits[roc])
                    grp.create_dataset('n_hits', data=self.NHits[roc])
                    grp.create_dataset('clusters', data=self.Clusters[roc])
                    grp.create_dataset('n_clusters', data=self.NClusters[roc])

    def add_data(self, data):
        info('adding data ... ')
        self.PBar.start(len(data))
        for i, event in enumerate(data):
            self.add_event(event)
            self.PBar.update(i)

    def add_event(self, event):
        if not event.pixels:
            return
        self.TriggerPhase.append(event.triggerPhases[0]) if event.triggerPhases else do_nothing()
        n_hits = zeros(self.NPlanes)
        for hit in event.pixels:
            n_hits[hit.roc] += 1
            self.Hits[hit.roc].append((hit.column, hit.row, hit.value, self.get_vcal(hit.roc, hit.column, hit.row, hit.value)))
        self.NHits.append(n_hits) if any(n_hits) else do_nothing()

    def make_arrays(self):
        for roc in range(self.NPlanes):
            self.Hits[roc] = array(self.Hits[roc], dtype=[('column', 'u2'), ('row', 'u2'), ('adc', 'i2'), ('vcal', 'f4')])
        self.NHits = array(self.NHits, 'u1').T
        self.NEvents = self.NHits[0].size

    # TODO revise this since algorithm is faulty!
    def clusterise(self):
        self.PBar.start(self.NEvents * self.NPlanes)
        info('clusterise ...')
        for roc in range(self.NPlanes):
            for i, event in enumerate(split(self.Hits[roc], cumsum(self.NHits[roc])[:-1])):
                cluster_hits = [[event[0]]]
                for hit in event[1:]:
                    belongs_to_existing_cluster = False
                    for clusters in cluster_hits:
                        # check if any hit is close to a cluster hit
                        if any((abs(array(clusters)['column'].astype('i2') - hit['column']) <= 1) & (abs(array(clusters)['row'].astype('i2') - hit['row']) <= 1)):
                            clusters.append(hit)
                            belongs_to_existing_cluster = True
                            break
                    # make new cluster
                    if not belongs_to_existing_cluster:
                        cluster_hits.append([hit])
                for cluster in cluster_hits:
                    cluster = array(cluster)
                    self.Clusters[roc].append((mean(cluster['column']), mean(cluster['row']), sum(cluster['vcal'])))
                self.NClusters[roc].append(len(cluster_hits))
                self.PBar.update(i + self.NEvents * roc)
            self.Clusters[roc] = array(self.Clusters[roc], dtype=[('column', 'f2'), ('row', 'f2'), ('vcal', 'f4')])
            self.NClusters[roc] = array(self.NClusters[roc], 'u1')


if __name__ == '__main__':

    z = HDF5Writer('main')
