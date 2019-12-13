#!/usr/bin/env python
# --------------------------------------------------------
#       Class to read pXar events from hdf5 file
# created on November 13th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

import h5py
from file_reader import *
from numpy import full, arange, where
from ROOT import TH2I, TH1F, TProfile2D, TProfile
from json import loads


class HDF5Reader(FileReader):

    def __init__(self, run_number=None, dut=0, config_name='main'):
        FileReader.__init__(self, run_number, config_name, file_type='hdf5')

        # Data
        self.Data = self.File['ROC{}'.format(dut)]
        self.NEvents = self.Data['n_hits'].size
        self.NHits = self.Data['hits'].size
        self.NClusters = self.Data['clusters'].size
        self.Fid = self.load_fiducial()
        self.FidCut = self.load_fid_cut()

        self.Bins = [self.NCols, arange(-.5, self.NCols), self.NRows, arange(-.5, self.NRows)]

    def __del__(self):
        self.File.close()
        
    def load_fid_cut(self):
        x, y = self.Data['clusters']['column'], self.Data['clusters']['row']
        return where((x >= self.Fid[0]) & (x <= self.Fid[1]) & (y >= self.Fid[2]) & (y <= self.Fid[3]))[0] if self.Fid else None

    def get_vcal_values(self, use_fid=True):
        return self.Data['clusters']['vcal'] if not use_fid or not self.Fid else self.Data['clusters']['vcal'][self.FidCut]
    
    def get_x(self, use_fid=True):
        return self.Data['clusters']['column'] if not use_fid or not self.Fid else self.Data['clusters']['column'][self.FidCut]

    def get_y(self, use_fid=True):
        return self.Data['clusters']['row'] if not use_fid or not self.Fid else self.Data['clusters']['row'][self.FidCut]

    def load_file(self):
        return h5py.File(join(self.DataDir, self.FileName))

    def draw_hitmap(self, cluster=False, vcal=None, fid=False):
        h = TH2I('hhm{}'.format(cluster), '{} Map'.format('Cluster' if cluster else 'Hit'), *self.Bins)
        if cluster:
            x = self.get_x(fid) if vcal is None else self.get_x(fid)[where(self.get_vcal_values(fid) < vcal)]
            y = self.get_y(fid) if vcal is None else self.get_y(fid)[where(self.get_vcal_values(fid) < vcal)]
            h.FillN(x.size, x.astype('d'), y.astype('d'), full(x.size, 1, 'd'))
        else:
            h.FillN(self.NClusters, self.Data['clusters']['column'].astype('d'), self.Data['clusters']['row'].astype('d'), full(self.NHits, 1, 'd'))
        format_histo(h, x_tit='column', y_tit='row', y_off=1.2, z_tit='Number of Entries', z_off=1.6, stats=0)
        self.Plotter.draw_histo(h, lm=.13, rm=.18, draw_opt='colz', x=1.17)

    def draw_cluster_map(self, vcal=None, fid=False):
        self.draw_hitmap(cluster=True, vcal=vcal, fid=fid)

    def draw_signal_map(self):
        h = TProfile2D('psm', 'Signal Map', *self.Bins)
        for x, y, v in zip(self.Data['clusters']['column'], self.Data['clusters']['row'], self.Data['clusters']['vcal']):
            h.Fill(x, y, v)
        format_histo(h, x_tit='column', y_tit='row', y_off=1.2, z_tit='VCAL', z_off=1.6, stats=0)
        self.Plotter.draw_histo(h, lm=.13, rm=.18, draw_opt='colz', x=1.17)

    @staticmethod
    def get_vcal_bins(bin_width):
        bins = arange(-200, 1800.01, bin_width)
        return [bins.size - 1, bins]

    def get_event_bins(self, bin_width=1000):
        bins = arange(0, (self.FidCut.size if self.Fid else self.NClusters) + .01, bin_width)
        return [bins.size - 1, bins]

    def draw_vcal(self, bin_width=5, use_fid=True):
        h = TH1F('hv', 'VCAL Distribution', *self.get_vcal_bins(bin_width))
        vcals = self.get_vcal_values(use_fid)
        h.FillN(vcals.size, vcals.astype('d'), full(vcals.size, 1, 'd'))
        v_max = h.GetMaximum()
        x_range = [h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(.01 * v_max) - 2, h.FindLastBinAbove(.01 * v_max) + 2]]
        format_histo(h, x_tit='VCAL', y_tit='Number of Entries', y_off=2, fill_color=self.Plotter.FillColor, x_range=x_range)
        self.Plotter.format_statbox(entries=True, m=True)
        self.Plotter.draw_histo(h, lm=.15)

    def draw_vcal_time(self, bin_width=1000):
        vcals = self.get_vcal_values()
        h = TProfile('hvt', 'VCAL vs. Time', *self.get_event_bins(bin_width))
        h.FillN(vcals.size, arange(vcals.size, dtype='d'), vcals.astype('d'), full(vcals.size, 1, 'd'))
        format_histo(h, x_tit='Cluster Number', y_tit='VCAL', y_off=1.5)
        self.Plotter.format_statbox(entries=1, x=.9)
        self.Plotter.draw_histo(h, lm=.12, rm=.08)

    def draw_trigger_phase(self):
        h = TH1F('htp', 'Trigger Phase Distribution', 10, 0, 10)
        h.FillN(self.NEvents, array(self.File['trigger_phase'], 'd'), full(self.NEvents, 1, 'd'))
        format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.2, fill_color=self.Plotter.FillColor)
        self.Plotter.format_statbox(entries=True)
        self.Plotter.draw_histo(h, lm=.12)

    def get_vcal(self):
        return mean_sigma(self.Data['clusters']['vcal'])

    def below_thresh(self, thresh=35):
        return sum(self.Data['clusters']['vcal'] < thresh) / float(self.NClusters) * 100

    def load_fiducial(self):
        if 'fid' in self.Config.options('CHIP'):
            return loads(self.Config.get('CHIP', 'fid'))


if __name__ == '__main__':

    from argparse import ArgumentParser
    aparser = ArgumentParser()
    aparser.add_argument('run', nargs='?', default=None, help='run number [default = the last saved run]')
    aparser.add_argument('-c', '--config', nargs='?', default='main', help='name of the ini config file (without .ini) [default = main]')
    pargs = aparser.parse_args()

    z = HDF5Reader(pargs.run, config_name=pargs.config)
