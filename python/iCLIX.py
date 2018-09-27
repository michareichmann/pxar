#!/usr/bin/env python
# --------------------------------------------------------
#       ipython command line tool using the pXar core api
# created on February 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from os.path import basename, dirname, realpath, join
from sys import argv, path, stdout

# add python and cython libs
pxar_dir = dirname(dirname(realpath(__file__)))
path.insert(1, join(pxar_dir, 'lib'))
path.insert(1, join(pxar_dir, 'python', 'src'))

from ROOT import TCanvas, TCutG, gStyle, TColor, TMultiGraph, TH1I
from argparse import ArgumentParser
from numpy import zeros, array, mean
from time import time, sleep
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from pxar_helpers import *
from pxar_plotter import Plotter
from TreeWriterLjubljana import TreeWriterLjubljana
from utils import *
from json import dumps

dacdict = PyRegisterDictionary()
probedict = PyProbeDictionary()
prog_name = basename(argv.pop(0))


class CLIX:
    """Simple command processor for the pxar core API."""

    def __init__(self, conf_dir, verbosity, trim):
        # main
        self.api = PxarStartup(conf_dir, verbosity, trim)
        self.IsAnalogue = 'dig' not in self.api.getRocType()
        self.Dir = conf_dir
        self.Trim = trim
        self.Verbosity = verbosity

        # dicts
        self.DacDict = PyRegisterDictionary()
        self.ProbeDict = PyProbeDictionary()
        self.TBDelays = self.api.getTestboardDelays()

        self.window = None
        self.Plots = []
        self.ProgressBar = None
        self.NRows = 80
        self.NCols = 52
        set_palette()

        self.Plotter = Plotter()

    def restart_api(self):
        self.api = None
        self.api = PxarStartup(self.Dir, self.Verbosity, self.Trim)

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()], maxval=n)
        self.ProgressBar.start()

    @staticmethod
    def run(filename):
        print '\nReading commands from file...\n'
        if not file_exists(filename):
            'File {} does not exit'.format(filename)
        f = open(filename)
        for line in f.readlines():
            if line.startswith('#'):
                continue
            print line.strip('\n\r')
            data = line.split()
            arguments = [float(word) if is_num(word) else word for word in data[1:]]
            exec 'z.{f}({args})'.format(f=data[0], args=dumps(arguments).strip('[]'))

    @staticmethod
    def remove_tbm_info(event):
        """Removes the TBM information (first 4bit) from the 16bit words."""
        return [word & 0x0fff for word in event]

    @staticmethod
    def expand_sign(event):
        """Retrieves the sign information of the 16bit words. ADC has only positive values. Value is negative if the third hex > 8"""
        return [word - (0x1000 if word & 0x0800 else 0) for word in event]

    def convert_raw_event(self, event):
        event = self.remove_tbm_info(event)
        return self.expand_sign(event) if self.IsAnalogue else event

    @staticmethod
    def decode_header(num):
        bin_str = bin(num)[2:]
        print 'Decoding Header:\n    MMMM 0111 1111 10RB'
        print 'bin {w}'.format(w=' '.join([bin_str[i:i + 4] for i in xrange(0, len(bin_str), 4)]))
        print 'hex    {w}'.format(w='    '.join(hex(num)[2:]))
        print 'header identifier: {hi} {eq} 0x7f8'.format(hi=hex(num & 0x0ffc), eq='=' if (num & 0x0ffc) == 0x7f8 else '!=')
        return (num & 0x0ffc) == 0x7f8

    @staticmethod
    def decode_digital(value):
        # 0x0fff0fff -> 0xffffff
        print '\nDecoding digital hit:\n   C1  C0  R2  R1  R0'
        bin_str = bin(value)[2:].zfill(6 * 4)
        print '   ' + ' '.join([bin_str[i:i + 3] for i in xrange(0, 5 * 3, 3)])
        print '   ' + ' '.join([' {} '.format(int(bin_str[i:i+3], 2)) for i in xrange(0, 5 * 3, 3)])
        print '\ncol = 2 * (6 * C1 + C0) + (R1 & 1)'
        print 'row = 80 - (36 * R2 + 6 * R1 + R0) / 2'
        col = 2 * (6 * bit_shift(value, 21) + bit_shift(value, 18)) + (bit_shift(value, 9) & 1)
        row = 80 - (36 * bit_shift(value, 15) + 6 * bit_shift(value, 12) + bit_shift(value, 9)) / 2
        ph = (value & 0x000f) + ((value >> 1) & 0x00f0)
        return row, col, ph

    def decode_pixel(self, lst):
        n_hits = 0
        for i in xrange(0, len(lst), 2):
            if (lst[i] & 0x0ffc) == 0x7f8:
                break
            print '\nDecoding Pixel Hit {n}'.format(n=i / 2 + 1)
            bin_str = ''.join(bin(lst[j])[2:].zfill(16) for j in [i, i + 1])
            print '    0000 CCCC CCRR RRRR MMMM RRRP PPP0 PPPP'
            print 'bin {w}'.format(w=' '.join([bin_str[j:j + 4] for j in xrange(0, len(bin_str), 4)]))
            print 'hex    {w}'.format(w='    '.join(hex(int(bin_str, 2))[2:].zfill(8)))
            raw_int = (lst[i] << 12) + (lst[i + 1] & 0x0fff)
            row, col, ph = self.decode_digital(raw_int) if not self.IsAnalogue else None
            print '\n===== [{c}, {r}, {p}] ====='.format(c=col, r=row, p=ph)
            n_hits += 2
            if hex(lst[i + 1])[:2].zfill(4).startswith('4'):
                break
        return n_hits

    # -----------------------------------------
    # region API
    def get_dac(self, dac, roc_id=0):
        dacs = self.api.getRocDACs(roc_id)
        if dac in dacs:
            return dacs[dac]
        else:
            print 'Unknown dac {d}!'.format(d=dac)

    def set_tb_delay(self, delay, value):
        if delay not in self.DacDict.getAllDTBNames():
            print 'The delay {} does not exist'.format(delay)
            return
        self.TBDelays[delay] = value
        self.api.setTestboardDelays(self.TBDelays)

    def set_dac(self, name, value, roc_id=None):
        self.api.setDAC(name, value, roc_id)

    def get_ia(self):
        print 'Analog Current: {} mA'.format(self.api.getTBia() * 1000)
    # endregion

    # -----------------------------------------
    # region DAQ
    def daq_start(self, arg=0):
        self.api.daqStart(arg)

    def daq_stop(self):
        self.api.daqStop()

    def daq_trigger(self, n_trig=1, period=500):
        self.api.daqTrigger(n_trig, period)

    def daq_get_event(self):
        try:
            data = self.api.daqGetEvent()
            print data
        except RuntimeError:
            pass

    def daq_get_raw_event(self, convert=1):
        try:
            event = self.api.daqGetRawEvent()
        except RuntimeError:
            return
        event = self.convert_raw_event(event) if convert == 1 else event
        return event

    def get_tb_ia(self):
        """ returns analog DTB current """
        print 'Analog Current: {c} mA'.format(c=self.api.getTBia() * 1000)

    def show_event_decoding(self):
        self.daq_start()
        self.daq_trigger()
        event = self.daq_get_raw_event(convert=0)
        print 'Raw data: {}\nhex:   {}\n'.format(event, [hex(num)[2:].zfill(4) for num in event])
        i = 0
        for _ in event:
            if i < len(event) and self.decode_header(event[i]):
                i += self.decode_pixel(event[i+1:])
            i += 1

    def signal_probe(self, probe=None, signal=None):
        probes = ['a1', 'a2', 'd1', 'd2']
        probe = raw_input('Enter the probe output {}: '.format(probes)) if probe is None else probe
        signals = self.ProbeDict.getAllAnalogNames() if probe.startswith('a') else self.ProbeDict.getAllDigitalNames()
        signal = raw_input('Enter a signal from {}: '.format(signals)) if signal is None else signal
        if probe not in probes or signal not in signals:
            print 'wrong probe or signal'
            return
        self.api.SignalProbe(probe, signal)

    def set_clock(self, value):
        """sets all the delays to the right value if you want to change clk"""
        self.set_tb_delay('clk', value)
        self.set_tb_delay('ctr', value)
        self.set_tb_delay('sda', value + (15 if 'dig' in self.api.getRocType() else 11))
        self.set_tb_delay('tin', value + (5 if 'dig' in self.api.getRocType() else 2))

    def set_pg(self, cal=True, res=True, trg=True, delay=None):
        """ Sets up the trigger pattern generator for ROC testing """
        pgcal = self.get_dac('wbc') + (6 if 'dig' in self.api.getRocType() else 5)
        pg_setup = []
        if delay is not None:
            pg_setup.append(('DELAY', delay))
        if res:
            pg_setup.append(('PG_RESR', 25))
        if cal:
            pg_setup.append(('PG_CAL', pgcal))
        if trg:
            pg_setup.append(('PG_TRG', 0 if self.api.getNTbms() != 0 else 15))
        if self.api.getNTbms() == 0:
            pg_setup.append(('PG_TOK', 0))
        # print pg_setup
        try:
            self.api.setPatternGenerator(tuple(pg_setup))
        except RuntimeError, err:
            print err

    def trigger_source(self, source, freq=0):
        """daqTriggerSource: select the trigger source to be used for the DAQ session"""
        if self.api.daqTriggerSource(source, 40000000 / freq if freq else 0):
            print 'Trigger source {} selected.'.format(source)
        else:
            print 'DAQ returns faulty state.'

    def trigger_loop(self, on='True', freq=100):
        """start\stop trigger loop: [on] [frequency]"""
        on = False if str(on).lower() in ['0', 'false', 'stop', 'end'] else True
        self.api.daqTriggerSource('periodic' if on else 'pg_dir', 40000000 / float(freq) if on else 0)
        self.daq_start()
        self.daq_trigger()
        self.daq_stop()
        print 'Trigger loop with frequency of {f}Hz {m}'.format(f=freq, m='started' if on else 'stopped')

    # endregion

    # -----------------------------------------
    # region MASK // ENABLE
    def get_activated(self, roc=None):
        return self.api.getNEnabledPixels(roc), self.api.getNMaskedPixels(roc)

    def enable_single_pixel(self, row=14, column=14, roc=None):
        """enableOnePixel [row] [column] [roc] : enables one Pixel (default 14/14); masks and disables the rest"""
        print '--> disable and mask all pixels of all activated ROCs'
        self.api.testAllPixels(0)
        self.api.maskAllPixels(1)
        self.api.testPixel(row, column, 1, roc)
        self.api.maskPixel(row, column, 0, roc)
        print_string = '--> enable and unmask Pixel {r}/{c}: '.format(r=row, c=column)
        print_string += '(' + ','.join('ROC {n}: {a}/{m}'.format(n=roc, a=self.get_activated(roc)[0], m=self.get_activated(roc)[1]) for roc in xrange(self.api.getNEnabledRocs())) + ')'
        print print_string

    def enable_all(self, roc=None):
        """enableAllPixel [roc]: enables and unmasks all Pixels of [roc]"""
        self.api.maskAllPixels(0, roc)
        self.api.testAllPixels(1, roc)
    # endregion

    # -----------------------------------------
    # region plotting
    def print_eff(self, data, n_trig):
        unmasked = 4160 - self.api.getNMaskedPixels()
        active = self.api.getNEnabledPixels()
        read_back = sum(px.value for px in data)
        total = n_trig * (unmasked if unmasked < active else active)
        eff = 100. * read_back / total
        print 'Efficiency: {eff:6.2f}% ({rb:5d}/{tot:5d})'.format(eff=eff, rb=int(read_back), tot=total)
        return eff

    def plot_graph(self, gr, lm=.12, rm=.1, draw_opt='alp'):
        c = TCanvas('c', 'c', 1000, 1000)
        c.SetMargin(lm, rm, .1, .1)
        gr.Draw(draw_opt)
        self.window = c
        self.Plots.append(gr)
        self.window.Update()

    def plot_map(self, data, name, count=False, no_stats=False):
        c = TCanvas('c', 'c', 1000, 1000)
        c.SetRightMargin(.12)

        # Find number of ROCs present:
        is_module = self.api.getNRocs() > 1
        proc = 'proc' in self.api.getRocType()
        # Prepare new numpy matrix:
        d = zeros((417 if is_module else 52, 161 if is_module else 80))
        for px in data:
            roc = (px.roc - 12) % 16 if proc else px.roc
            xoffset = 52 * (roc % 8) if is_module else 0
            yoffset = 80 * int(roc / 8) if is_module else 0

            # Flip the ROCs upside down:
            y = (px.row + yoffset) if (roc < 8) else (2 * yoffset - px.row - 1)
            # Reverse order of the upper ROC row:
            x = (px.column + xoffset) if (roc < 8) else (415 - xoffset - px.column)
            d[x][y] += 1 if count else px.value

        plot = Plotter.create_th2(d, 0, 417 if is_module else 52, 0, 161 if is_module else 80, name, 'pixels x', 'pixels y', name)
        if no_stats:
            plot.SetStats(0)
        plot.Draw('COLZ')
        # draw margins of the ROCs for the module
        self.draw_module_grid(is_module)
        self.window = c
        self.Plots.append(plot)
        self.window.Update()

    def draw_module_grid(self, draw):
        if not draw:
            return
        for i in xrange(2):
            for j in xrange(8):
                rows, cols = self.NRows, self.NCols
                x = array([cols * j, cols * (j + 1), cols * (j + 1), cols * j, cols * j], 'd')
                y = array([rows * i, rows * i, rows * (i + 1), rows * (i + 1), rows * i], 'd')
                cut = TCutG('r{n}'.format(n=j + (j * i)), 5, x, y)
                cut.SetLineColor(1)
                cut.SetLineWidth(1)
                self.Plots.append(cut)
                cut.Draw('same')
    # endregion

    def get_efficiency_map(self, flags=0, n_triggers=10):
        data = self.api.getEfficiencyMap(flags, n_triggers)
        self.print_eff(data, n_triggers)
        self.plot_map(data, 'Efficiency', no_stats=True)

    def clk_scan(self, exclude=None):
        """ scanning digital clk and deser phases """
        n = 10
        n_rocs = self.api.getNRocs()
        self.set_pg(cal=False, res=True)
        self.daq_start()
        print '\nCLK',
        for i in xrange(8):
            print '{:2d} '.format(i),
        print
        good = []
        for clk in xrange(20):
            if clk == exclude:
                continue
            self.set_clock(clk)
            print '{:2d}:'.format(clk),
            for phase in xrange(8):
                self.set_tb_delay('deser160phase', phase)
                # self.api.setTestboardDelays({'clk': clk})
                self.daq_trigger(n)
                evts = [self.daq_get_raw_event() for _ in xrange(n)]
                eff = mean([1 if event is not None and len(event) == n_rocs and all(header in xrange(2040, 2044) for header in event) else 0 for event in evts])
                if eff == 1:
                    good.append((clk, phase))
                    print '{c}{eff:1.1f}{e}'.format(eff=eff, c=GREEN, e=ENDC),
                elif eff > .5:
                    print '{c}{eff:1.1f}{e}'.format(eff=eff, c=YELLOW, e=ENDC),
                elif eff > 0:
                    print '{c}{eff:1.1f}{e}'.format(eff=eff, c=RED, e=ENDC),
                else:
                    print ' x ',
            print
        self.daq_stop()
        self.set_pg(cal=True, res=True)
        if not good:
            print 'Did not find any good timing...'
            return
        clk, phase = good[(len(good) / 2)]
        print 'Set CLK/DESER160PHASE to: {}/{}'.format(clk, phase)
        self.set_clock(clk)
        self.set_tb_delay('deser160phase', phase)

    def scan_clk(self):
        self.daq_start()
        for clk in xrange(20):
            self.daq_trigger()
            self.set_clock(clk)
            print '{:2d}:'.format(clk), self.daq_get_raw_event()
        self.daq_stop()

    def do_adc_disto(self, vcal=50, col=14, row=14, high=False, n_trig=10000):
        self.api.setDAC('ctrlreg', 4 if high else 0)
        self.api.setDAC('vcal', vcal)
        self.api.daqStart()
        self.api.daqTrigger(n_trig, 500)
        self.enable_single_pixel(col, row)
        data = self.api.daqGetEventBuffer()
        self.api.daqStop()
        adcs = [px.value for evt in data for px in evt.pixels]
        h = TH1I('h_adc', 'ACD Distribution for vcal {v} in {h} Range'.format(v=vcal, h='high' if high else 'low'), 255, 0, 255)
        for adc in adcs:
            h.Fill(adc)
        self.plot_graph(h, draw_opt='')

    def clear_buffer(self):
        try:
            self.api.daqGetEventBuffer()
        except RuntimeError:
            pass

    def wbc_scan(self, min_wbc=90, max_triggers=50, max_wbc=130):
        """do_wbcScan [minimal WBC] [number of events] [maximal WBC]: \n
        sets wbc from minWBC until it finds the wbc which has more than 90% filled events or it reaches maxWBC \n
        (default [90] [100] [130])"""

        # prepararations
        print 'Turning on HV!'
        self.api.HVon()
        print 'Setting trigger source to "extern"'
        self.api.daqTriggerSource('extern')
        self.api.daqStart()

        trigger_phases = zeros(10)
        yields = OrderedDict([(roc, {wbc: 0 for wbc in xrange(min_wbc, max_wbc)}) for roc in xrange(self.api.getNEnabledRocs())])
        set_dacs = {roc: False for roc in yields}
        print '\nROC EVENT YIELDS:\n  wbc\t{r}'.format(r='\t'.join(('roc' + str(yld)).rjust(6) for yld in yields.keys()))

        # loop over wbc
        for wbc in xrange(min_wbc, max_wbc):
            self.clear_buffer()
            self.api.setDAC('wbc', wbc)

            # loop until you find nTriggers
            n_triggers = 0
            while n_triggers < max_triggers:
                try:
                    data = self.api.daqGetEvent()
                    trigger_phases[data.triggerPhases[0]] += 1
                    for roc in set([pix.roc for pix in data.pixels]):
                        yields[roc][wbc] += 1. * 100. / max_triggers
                    n_triggers += 1
                except RuntimeError:
                    pass
            y_strings = ['{y:5.1f}%'.format(y=yld) for yld in [yields[roc][wbc] for roc in yields.iterkeys()]]
            print '  {w:03d}\t{y}'.format(w=wbc, y='\t'.join(y_strings))

            # stopping criterion
            best_wbc = max_wbc
            if wbc > min_wbc + 3:
                for roc, ylds in yields.iteritems():
                    if any(yld > 10 for yld in ylds.values()[:wbc - min_wbc - 2]):
                        best_wbc = ylds.keys()[ylds.values().index(max(ylds.values()))]
                        print 'set wbc of roc {i} to {v}'.format(i=roc, v=best_wbc)
                        self.api.setDAC('wbc', best_wbc, roc)
                        set_dacs[roc] = True
                if all(set_dacs.itervalues()):
                    print 'found all wbcs'

                    for roc, dic in yields.iteritems():
                        keys = dic.keys()
                        for key in keys:
                            if key >= best_wbc + 4:
                                yields[roc].pop(key)
                    break
        self.api.daqStop()

        # triggerphase

        print '\nTRIGGER PHASE:'
        for i, trigger_phase in enumerate(trigger_phases):
            if trigger_phase:
                percentage = trigger_phase * 100 / sum(trigger_phases)
                print '{i}\t{d} {v:2.1f}%'.format(i=i, d=int(round(percentage)) * '|', v=percentage)

        # plot wbc_scan
        mg = TMultiGraph('mg_wbc', 'WBC Scans for all ROCs')
        colors = range(1, len(yields) + 1)
        leg = Plotter.create_legend(nentries=len(yields), x1=.7)
        for i, (roc, dic) in enumerate(yields.iteritems()):
            gr = Plotter.create_graph(x=dic.keys(), y=dic.values(), tit='wbcs for roc {r}'.format(r=roc), xtit='wbc', ytit='yield [%]', color=colors[i])
            leg.AddEntry(gr, 'roc{r}'.format(r=roc), 'lp')
            mg.Add(gr, 'lp')
        self.Plotter.plot_histo(mg, draw_opt='a', l=leg)
        mg.GetXaxis().SetTitle('WBC')
        mg.GetYaxis().SetTitle('Yield [%]')

    def hitmap(self, t=1, random_trigger=1, n=10000):
        self.api.HVon()
        t_start = time()
        if random_trigger:
            self.set_pg(cal=False, res=False, delay=20)
        else:
            self.signal_probe('a1', 'sdata2')
            self.set_dac('wbc', 119)
            self.api.daqTriggerSource('extern')
        self.api.daqStart()
        self.start_pbar(t * 600)
        data = []
        while time() - t_start < t * 60:
            self.ProgressBar.update(int((time() - t_start) * 10) + 1)
            if random_trigger:
                self.api.daqTrigger(n, 500)
            try:
                if random_trigger:
                    sleep(.5)
                    data += self.api.daqGetEventBuffer()
                else:
                    data += [self.api.daqGetEvent()]
            except RuntimeError:
                pass
            except MemoryError:
                break
        self.ProgressBar.finish()
        self.api.daqStop()
        self.api.HVoff()
        self.set_pg()
        pix_data = [pix for event in data for pix in event.pixels]
        h = TH1I('h', 'h', 512, -256, 256)
        for pix in pix_data:
            h.Fill(pix.value)
        print 'Entries:', h.GetEntries
        self.Plotter.plot_histo(h, draw_opt='hist')
        self.plot_map(pix_data, 'Hit Map', count=True, no_stats=True)
        stats = self.api.getStatistics()
        event_rate = stats.valid_events / (2.5e-8 * stats.total_events / 8.)
        hit_rate = stats.valid_pixels / (2.5e-8 * stats.total_events / 8.)
        stats.dump
        print 'Event Rate: {0:5.4f} MHz'.format(event_rate / 1000000)
        print 'Hit Rate:   {0:5.4f} MHz'.format(hit_rate / 1000000)

    def hitmap_random(self, t, n=10000):
        return self.hitmap(t, random_trigger=True, n=n)

    def hitmap_trigger(self, t):
        return self.hitmap(t, random_trigger=False)

    def load_mask(self, file_name):
        f = open(file_name, 'r')
        lines = [line.strip('\n') for line in f.readlines() if len(line) > 3 and not line.startswith('#')]
        for i, line in enumerate(lines):
            data1 = line.split(' ')
            if data1[0] == 'cornBot':
                data2 = lines[i + 1].split(' ')
                i2c = int(data1[1])
                self.api.maskAllPixels(True, i2c)
                for col in xrange(int(data1[2]), int(data2[2]) + 1):
                    for row in xrange(int(data1[3]), int(data2[3]) + 1):
                        print col, row, i2c
                        self.api.maskPixel(col, row, False, i2c)

    def save_data(self):
        t = TreeWriterLjubljana()
        self.api.HVon()
        self.trigger_source('extern')
        self.set_dac('wbc', t.Config.getint('MAIN', 'wbc'))
        self.signal_probe('a1', 'sdata2')
        self.daq_start()
        i = 0
        while True:
            try:
                t.write(self.api.daqGetEvent())
                print '\r{}'.format(i),
                stdout.flush()
                i += 1
            except RuntimeError:
                pass
            except KeyboardInterrupt:
                info('Data taking stopped by KeyboardInterrupt!')
                break
        self.daq_stop()


def set_palette(custom=True, pal=1):
    if custom:
        stops = array([0., .5, 1], 'd')
        green = array([0. / 255., 200. / 255., 80. / 255.], 'd')
        blue = array([0. / 255., 0. / 255., 0. / 255.], 'd')
        red = array([180. / 255., 200. / 255., 0. / 255.], 'd')
        gStyle.SetNumberContours(20)
        bla = TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 255)
        color_table = array([bla + ij for ij in xrange(255)], 'i')
        gStyle.SetPalette(len(color_table), color_table)
    else:
        gStyle.SetPalette(pal)


def do_nothing():
    pass


def bit_shift(value, shift):
    return (value >> shift) & 0b0111


if __name__ == '__main__':
    # command line argument parsing

    parser = ArgumentParser(prog=prog_name, description="A Simple Command Line Interface to the pxar API.")
    parser.add_argument('--dir', '-d', metavar="DIR", help="The digit rectory with all required config files.")
    parser.add_argument('--run', '-r', metavar="FILE", help="Load a cmdline script to be executed before entering the prompt.", default='')
    parser.add_argument('--verbosity', '-v', metavar="LEVEL", default="INFO", help="The output verbosity set in the pxar API.")
    parser.add_argument('--trim', '-T', nargs='?', default=None, help="The output verbosity set in the pxar API.")
    parser.add_argument('-wbc', action='store_true')
    args = parser.parse_args(argv)

    print_banner('# STARTING ipython pXar Command Line Interface')

    # start command line
    z = CLIX(args.dir, args.verbosity, args.trim)
    if args.wbc:
        z.wbc_scan()
        raw_input('Enter any key to close the program')
    print
    if args.run:
        z.run(args.run)

    # shortcuts

    ia = z.get_ia
    hvon = z.api.HVon
    hvoff = z.api.HVoff
    ge = z.get_efficiency_map
    ds = z.daq_start
    st = z.daq_stop
    ev = z.daq_get_event
    raw = z.daq_get_raw_event
    dt = z.daq_trigger
    # sd = z.api.setDAC