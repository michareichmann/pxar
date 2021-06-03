#!/usr/bin/env python
# --------------------------------------------------------
#       ipython command line tool using the pXar core api
# created on February 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from os.path import basename, dirname, realpath, join
from subprocess import call
from sys import argv, path, stdout
from os import _exit as terminate

# add python and cython libs
pxar_dir = dirname(dirname(realpath(__file__)))
path.insert(1, join(pxar_dir, 'lib'))
path.insert(1, join(pxar_dir, 'python', 'src'))

from argparse import ArgumentParser
from numpy import delete, argmin, sign, argmax, genfromtxt
from numpy.random import randint
from time import sleep
from pxar_helpers import *
from TreeWriterLjubljana import TreeWriterLjubljana
from hdf5_writer import HDF5Writer
from json import dumps
from threading import Thread
import atexit
from draw import *

BREAK = False
z = None
import signal
def signal_handler(signal, frame):
    global BREAK, z
    print("\nprogram exiting gracefully")
    BREAK = True
    terminate(5)


signal.signal(signal.SIGINT, signal_handler)
dacdict = PyRegisterDictionary()
probedict = PyProbeDictionary()
prog_name = basename(argv.pop(0))


def ex():
    if z is not None:
        z.IsRunning = False
        sleep(.5)
        z.api.HVoff()
    print 'Bye ...'


atexit.register(ex)


class CLIX:
    """Simple command processor for the pxar core API."""

    def __init__(self, conf_dir, verbosity, trim):
        # main
        self.api = PxarStartup(conf_dir, verbosity, trim)
        self.Config = PxarConfigFile(join(conf_dir, 'configParameters.dat'))
        self.NRocs, self.I2C = find_i2cs(self.Config, prnt=False)
        self.IsAnalogue = 'dig' not in self.api.getRocType()
        self.Dir = conf_dir
        self.Trim = trim
        self.Verbosity = verbosity

        # dicts
        self.DacDict = PyRegisterDictionary()
        self.ProbeDict = PyProbeDictionary()
        self.TBDelays = self.api.getTestboardDelays()

        self.NCols, self.NRows = 52, 80
        self.Draw = Draw()
        self.PBar = PBar()
        self.IsRunning = False

        self.get_ia()

    def restart_api(self):
        self.api = None
        self.api = PxarStartup(self.Dir, self.Verbosity, self.Trim)

    # -----------------------------------------
    # region HELPERS
    @staticmethod
    def run(filename):
        print '\nReading commands from file...\n'
        if not file_exists(filename):
            'File {} does not exist'.format(filename)
        f = open(filename)
        for line in f.readlines():
            if line.startswith('#'):
                continue
            data = line.split()
            arguments = [float(word) if is_num(word) else word for word in data[1:]]
            print 'z.{f}({args})'.format(f=data[0], args=dumps(arguments).strip('[]'))
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
    def translate_level(value, level1, level_s, offset=0):
        level = int((value + level1 + level_s + offset) / level1)
        return 5 if level >= 5 else level

    def decode_analogue(self, value, offset=0):
        black = value[1]
        level1 = (black - value[0]) / 4
        level_s = level1 / 2
        ph = value[-1]

        c1 = self.translate_level(value[3], level1, level_s, offset)
        c0 = self.translate_level(value[4], level1, level_s, offset)

        r2 = self.translate_level(value[5], level1, level_s, offset)
        r1 = self.translate_level(value[6], level1, level_s, offset)
        r0 = self.translate_level(value[7], level1, level_s, offset)

        column, row = calculate_col_row(c1, c0, r2, r1, r0)

        print ' S,  c1  c0  r2  r1  r0'
        print '{:2d}, {}'.format(level_s, ' '.join(['{:3d}'.format(i) for i in value[3:8]]))
        print '{:2d}, {}\n'.format(level_s, ' '.join(['{:3d}'.format(i) for i in [c1, c0, r2, r1, r0]]))
        print '\n===== [{c}, {r}, {p}] ====='.format(c=column, r=row, p=ph)

    def set_offset(self, value):
        self.api.setBlackOffsets(make_list(value))
        print 'set analogue decoding offset to: {}'.format(value)

    def set_l1_offset(self, value):
        self.api.setDecodingL1Offsets(make_list(value))
        print 'set analogue L1 decoding offset to: {}'.format(value)

    def set_alphas(self, value):
        self.api.setDecodingAlphas(make_list(value))
        print 'set analogue decoding alphas to: {}'.format(value)

    @staticmethod
    def decode_digital(value):
        # 0x0fff0fff -> 0xffffff
        print '\nDecoding digital hit:\n   C1  C0  R2  R1  R0'
        bin_str = bin(value)[2:].zfill(6 * 4)
        print '   ' + ' '.join([bin_str[i:i + 3] for i in xrange(0, 5 * 3, 3)])
        print '   ' + ' '.join([' {} '.format(int(bin_str[i:i+3], 2)) for i in xrange(0, 5 * 3, 3)])
        print '\ncol = 2 * (6 * C1 + C0) + (R0 & 1)'
        print 'row = 80 - (36 * R2 + 6 * R1 + R0) / 2'
        col, row = calculate_col_row(bit_shift(value, 21), bit_shift(value, 18), bit_shift(value, 15), bit_shift(value, 12), bit_shift(value, 9))
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

    @staticmethod
    def save_parameters(filename, current):
        data = genfromtxt(filename, dtype=['i', 'U30', 'i'])
        for i, (_, key, value) in enumerate(data):
            if key in current and current[key] != value:
                info('overwriting {} from {} -> {}'.format(key, value, current[key]))
                data[i][-1] = current[key]
        s = [max(len(str(i)) for i in data[a]) for a in ['f0', 'f1', 'f2']]
        with open(filename, 'r+') as f:
            f.writelines(['{0:>{3}}  {1:>{4}}  {2:>{5}}\n'.format(*list(tup) + s) for tup in data])
            f.truncate()

    def save_tb_delays(self):
        return self.save_parameters(self.Config.get('tbParameters'), self.TBDelays)

    def save_dac_parameters(self, roc_id=0, trim=None):
        self.save_parameters(self.get_dacfile(roc_id, trim), self.api.getRocDACs(roc_id))

    def save_config(self, key, value, name='configParameters'):
        f_name = join(self.Dir, '{}.dat'.format(name.split('.')[0]))
        with open(f_name, 'r+') as f:
            value = str(value)
            lines = f.readlines()
            for i, line in enumerate(lines):
                if key in line:
                    old = line[line.find(key) + len(key):]
                    rsize = len(old) - 1 if abs(len(value) - len(old.strip())) < 3 else len(value) + 1  # -1 for \n
                    lines[i] = line.replace(old, value.rjust(rsize)) + '\n'
                    info('changing existing config setting of {} from {} to {}'.format(key, ' '.join(line.split()[1:]), value))
            if not any(key in line for line in lines):
                lines.append('{} {}\n'.format(key, value))
                info('adding config setting {} with value {}'.format(key, value))
            f.seek(0)
            f.truncate()
            f.writelines(lines)

    def get_dacfile(self, roc_id=0, trim=None):
        return '{}{}_C{}.dat'.format(self.Config.get('dacParameters'), choose(choose(trim, self.Trim), ''), self.I2C[roc_id])
    # endregion HELPERS
    # -----------------------------------------

    # -----------------------------------------
    # region API
    def get_dac(self, dac, roc_id=0):
        """:returns: the current value of the DAC [dac] of the ROC with id [roc_id]. """
        dacs = self.api.getRocDACs(roc_id)
        if dac in dacs:
            return dacs[dac]
        else:
            print 'Unknown dac {d}!'.format(d=dac)

    def set_dac(self, name, value, roc_id=None):
        """sets the value of the DAC [dac] with ROC ID [roc_id] to [value]. """
        self.api.setDAC(name, value, roc_id)

    def get_tb_delay(self, name):
        """:returns: the current value of the testboard delay [name]. """
        delays = self.api.getTestboardDelays()
        if name not in self.api.getTestboardDelays():
            warning('{} not a valid delay'.format(name))
            return
        return delays[name]

    def set_tb_delay(self, delay, value):
        """sets the value of the DAC [dac] to [value]. """
        if delay not in self.DacDict.getAllDTBNames():
            print 'The delay {} does not exist'.format(delay)
            return
        self.TBDelays[delay] = value
        self.api.setTestboardDelays(self.TBDelays)

    def set_tb_delays(self, dic):
        for key, value in dic.items():
            self.set_tb_delay(key, value)

    def read_ia(self, t=.05):
        sleep(t)
        return self.api.getTBia()

    def get_ia(self, n=5, prnt=True):
        """:returns: the analogue current consumption of the testboard."""
        self.api.getTBia()  # first reading is always wrong
        current = mean([self.read_ia() for _ in range(n)]) * 1000  # to mA
        info('analogue current: {:.2f} mA'.format(current), prnt=prnt)
        return current if not prnt else None

    def get_n_rocs(self):
        """:returns: the number of enabled ROCs."""
        return self.api.getNEnabledRocs()
    # endregion API
    # -----------------------------------------

    # -----------------------------------------
    # region DAQ
    def clear_buffer(self):
        self.api.daqClear()

    def daq_start(self, arg=0):
        """starts the data acquisition."""
        self.api.daqStart(arg)

    def daq_stop(self):
        self.api.daqStop()

    def cycle_tb(self):
        self.daq_start()
        self.daq_stop()

    def daq_trigger(self, n_trig=1, period=500):
        self.api.daqTrigger(n_trig, period)

    def get_event(self):
        try:
            return self.api.daqGetEvent()
        except RuntimeError:
            return

    def update_time(self, t_start, t_max, n, n_max):
        t = time() - t_start
        if t - self.PBar.LastUpdate > .1:
            self.PBar.update(int(t * 10) if n_max is None else n)
            self.PBar.LastUpdate = t
        return t < t_max * 60 if n_max is None else n < n_max

    def take_data(self, wbc, t, n=None, random_trig=False):
        self.api.HVon()
        if random_trig:
            self.set_pg(cal=False, res=False, delay=20)
        else:
            self.signal_probe('a1', 'sdata2')
            self.set_dac('wbc', wbc)
            self.api.daqTriggerSource('extern')
        self.PBar.start(t * 600 if n is None else n)
        data = []
        self.daq_start()
        t_start = time()
        while time() - t_start < t * 60 if n is None else len(data) < n:
            sleep(.01)
            self.daq_trigger(10000) if random_trig else do_nothing()
            try:
                self.set_dac('wbc', wbc)  # resets the ROC ... lazy solution
                while self.update_time(t_start, t, len(data), n):
                    data.append(self.api.daqGetEvent())
            except RuntimeError:
                pass
        self.PBar.finish()
        self.print_rate(time() - t_start, random_trig)
        self.daq_stop()
        self.api.HVoff()
        self.set_pg()
        return data

    def print_rate(self, t, random_trig):
        stats = self.api.getStatistics()
        t = (2.5e-8 * stats.total_events / 8.) if random_trig else t
        stats.dump()
        print 'Event Rate: {0: 8.4f} kHz'.format(stats.valid_events / t / 1e3)
        print 'Hit Rate:   {0: 8.4f} kHz'.format(stats.valid_pixels / t / 1e3)
        print 'Trigger Eff {0: 8.4f} %'.format(100. * stats.valid_events / float(stats.total_events))

    def get_event_data(self, n):
        data = []
        while len(data) < n:
            event = self.get_event()
            data.append(event) if event is not None else do_nothing()
        return data

    def daq_get_raw_event(self, convert=True):
        try:
            event = self.api.daqGetRawEvent()
        except RuntimeError:
            return
        event = self.convert_raw_event(event) if convert else event
        return event

    def send_triggers(self, n, period=500):
        self.daq_start()
        self.daq_trigger(n, period)
        sleep(.2)
        self.daq_stop()
        sleep(.2)

    def get_raw_buffer(self, convert=True):
        events = []
        while not events or events[-1] is not None:
            events.append(self.daq_get_raw_event(convert))
        return array(events[:-1])

    def get_raw_event(self, convert=True, trigger=True, n_trig=1):
        if trigger:
            self.send_triggers(n_trig)
        return self.get_raw_buffer(convert)

    def get_mean_black(self, n_trigger=1000):
        self.send_triggers(n_trigger)
        return mean_sigma(self.get_raw_buffer()[:, 1])

    def get_tb_ia(self):
        """ returns analog DTB current """
        print 'Analog Current: {c} mA'.format(c=self.api.getTBia() * 1000)

    def show_event_decoding(self):
        self.daq_start()
        self.daq_trigger()
        event = self.daq_get_raw_event(convert=False)
        print 'Raw data: {}\nhex:   {}\n'.format(event, [hex(num)[2:].zfill(4) for num in event])
        i = 0
        for _ in event:
            if i < len(event) and self.decode_header(event[i]):
                i += self.decode_pixel(event[i+1:])
            i += 1

    def signal_probe(self, probe=None, sig=None):
        probes = ['a1', 'a2', 'd1', 'd2']
        probe = raw_input('Enter the probe output {}: '.format(probes)) if probe is None else probe
        signals = self.ProbeDict.getAllAnalogNames() if probe.startswith('a') else self.ProbeDict.getAllDigitalNames()
        sig = raw_input('Enter a signal from {}: '.format(signals)) if sig is None else sig
        if probe not in probes or sig not in signals:
            print 'wrong probe or signal'
            return
        return self.api.SignalProbe(probe, sig)

    @update_pbar
    def set_clock(self, value):
        """sets all the delays to the right value if you want to change clk"""
        self.set_tb_delay('clk', value)
        self.set_tb_delay('ctr', value)
        self.set_tb_delay('sda', value + (15 if 'dig' in self.api.getRocType() else 11))
        self.set_tb_delay('tin', value + (5 if 'dig' in self.api.getRocType() else 2))

    def set_external_clock(self, status=True):
        """setExternalClock [status]: enables the external DTB clock input, switches off the internal clock. Only switches if external clock is present."""
        print 'using {} clock {}'.format('external' if status else 'internal', 'failed' if not self.api.setExternalClock(status) else '')

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
        """start/stop trigger loop: [on] [frequency]"""
        on = False if str(on).lower() in ['0', 'false', 'stop', 'end'] else True
        self.api.daqTriggerSource('periodic' if on else 'pg_dir', 40000000 / float(freq) if on else 0)
        self.daq_start()
        self.daq_trigger()
        self.daq_stop()
        print 'Trigger loop with frequency of {f}Hz {m}'.format(f=freq, m='started' if on else 'stopped')

    def data_loop(self, freq=10, status=ON):

        if self.IsRunning:
            self.IsRunning = False
            sleep(.1)
        self.IsRunning = status
        if status == OFF or not freq:
            self.daq_stop()
            return

        def dloop():
            t0 = time()
            self.daq_start()
            z.set_clock(z.get_tb_delay('clk'))  # reduces noise
            while self.IsRunning:
                self.daq_trigger(1, 500)
                while time() - t0 < .99 / freq:
                    sleep(1. / 100. / freq)
                t0 = time()
        t = Thread(target=dloop)
        t.setDaemon(True)
        t.start()

    def get_raw_data(self, n_trigger=1000):
        self.send_triggers(n_trigger)
        return self.get_raw_buffer()

    def get_data(self, n_trigger=1000):
        self.send_triggers(n_trigger)
        return self.api.daqGetEventBuffer()

    def get_address_levels(self, n_trigger=1000):
        data = self.get_raw_data(n_trigger)[:, 3:]  # remove the ROC header
        return delete(data, range(5, data.shape[1], 6), axis=1)  # remove every sixth column (ph)

    def get_header(self, n_trigger=1000):
        data = self.get_raw_data(n_trigger)
        return data[:, argmin(data[0]):argmin(data[0]) + 3] if len(data.shape) == 2 else zeros((1, 3))  # take only header

    def get_mean_address_levels(self, n_trigger=1000):
        return mean(self.get_address_levels(n_trigger), axis=0)

    def get_mean_header(self, clk=None, n_trigger=1000):
        if clk is not None:
            self.set_clock(clk)
        return mean(self.get_header(n_trigger), axis=0)
    # endregion DAQ
    # -----------------------------------------

    # -----------------------------------------
    # region MASK // ENABLE
    def get_activated(self, roc=None):
        return self.api.getNEnabledPixels(roc), self.api.getNMaskedPixels(roc)

    def disable_all(self, roc=None):
        self.api.testAllPixels(0, roc)
        self.api.maskAllPixels(1, roc)

    def enable_pix(self, column, row, roc=None):
        self.api.testPixel(column, row, 1, roc)
        self.api.maskPixel(column, row, 0, roc)

    def enable_pixels(self, columns, rows, roc=None):
        for column, row in zip(columns, rows):
            self.enable_pix(column, row, roc)

    def enable_single_pixel(self, column=14, row=14, roc=None, prnt=True):
        """enableOnePixel [row] [column] [roc] : enables one Pixel (default 14/14); masks and disables the rest"""
        self.disable_all(roc)
        self.enable_pix(column, row, roc)
        if prnt:
            print '--> disable and mask all pixels of all activated ROCs'
            print_string = '--> enable and unmask Pixel {c}/{r}: '.format(r=row, c=column)
            print_string += '(' + ','.join('ROC {n}: {a}/{m}'.format(n=roc, a=self.get_activated(roc)[0], m=self.get_activated(roc)[1]) for roc in xrange(self.api.getNEnabledRocs())) + ')'
            print print_string

    def enable_all(self, roc=None):
        """enableAllPixel [roc]: enables and unmasks all Pixels of [roc]"""
        self.api.maskAllPixels(0, roc)
        self.api.testAllPixels(1, roc)
    # endregion MASK // ENABLE
    # -----------------------------------------

    # -----------------------------------------
    # region PLOTTING
    def print_eff(self, data, n_trig):
        unmasked = 4160 * self.get_n_rocs() - self.api.getNMaskedPixels()
        active = self.api.getNEnabledPixels()
        read_back = sum(px.value for px in data)
        total = n_trig * (unmasked if unmasked < active else active)
        eff = 100. * read_back / total
        print 'Efficiency: {eff:6.2f}% ({rb:5d}/{tot:5d})'.format(eff=eff, rb=int(read_back), tot=total)
        return eff

    def plot_map(self, data, title, count=False, stats=True):
        is_module = self.NRocs > 1
        proc = 'proc' in self.api.getRocType()
        x, y, zz = [], [], []  # Prepare new numpy matrix:
        for px in data:
            roc = (px.roc - 12) % 16 if proc else px.roc
            xoffset = 52 * (roc % 8) if is_module else 0
            yoffset = 80 * int(roc / 8) if is_module else 0
            y.append(px.row + yoffset if roc < 8 else 2 * yoffset - px.row - 1)  # Flip the ROCs upside down:
            x.append(px.column + xoffset if roc < 8 else 415 - xoffset - px.column)  # Reverse order of the upper ROC row:
            zz.append(px.value)
        if not count:
            x, y = [array(arr).repeat(zz) for arr in [x, y]]
        binning = make_bins(0, 417 if is_module else 52) + make_bins(0, 161 if is_module else 80)
        self.Draw.histo_2d(x, y, binning, title, x_tit='col', y_tit='row', stats=stats, z_range=[0, max(zz)])
        self.draw_module_grid(is_module)

    def draw_module_grid(self, draw):
        if draw:
            return [self.Draw.box(self.NCols * j, self.NRows * i, self.NCols * (j + 1), self.NRows * (i + 1)) for i in range(2) for j in range(8)]
    # endregion PLOTTING
    # -----------------------------------------

    def get_efficiency_map(self, flags=0, n_triggers=10):
        data = self.api.getEfficiencyMap(flags, n_triggers)
        self.print_eff(data, n_triggers)
        self.plot_map(data, 'Efficiency Map', stats=False)

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

    def draw_adc_disto(self, vcal=50, col=14, row=14, high=False, n_trig=10000):
        self.api.setDAC('ctrlreg', 4 if high else 0)
        self.api.setDAC('vcal', vcal)
        self.enable_single_pixel(col, row)
        x = [px.value for evt in self.get_data(n_trig) for px in evt.pixels]
        self.Draw.distribution(x, make_bins(-256, 256), 'ACD Distribution for vcal {} in {} Range'.format(vcal, 'high' if high else 'low'), x_tit='ADC', x_range=ax_range(x, 0, .2, .7))

    def wbc_scan(self, min_wbc=97, max_triggers=50, max_wbc=130, plot=False):
        """do_wbcScan [minimal WBC] [number of events] [maximal WBC]: \n
        sets wbc from minWBC until it finds the wbc which has more than 90% filled events or it reaches maxWBC \n
        (default [90] [100] [130])"""

        # preparations
        print 'Turning on HV!'
        self.api.HVon()
        print 'Setting trigger source to "extern"'
        self.api.daqTriggerSource('extern')
        self.api.SignalProbe('a1', 'sdata2')
        self.api.daqStart()

        trigger_phases = zeros(10)
        yields = OrderedDict((wbc, zeros(self.get_n_rocs())) for wbc in xrange(min_wbc, max_wbc))
        print '\nROC EVENT YIELDS:\n  wbc\t{r}'.format(r='\t'.join(('roc{}'.format(i)).rjust(6) for i in xrange(self.get_n_rocs())))

        for wbc in xrange(min_wbc, max_wbc):  # loop over wbc
            self.clear_buffer()
            self.api.setDAC('wbc', wbc)
            for event in self.get_event_data(max_triggers):
                if len(event.triggerPhases):
                    trigger_phases[event.triggerPhases[0]] += 1
                    for roc in set([pix.roc for pix in event.pixels]):
                        yields[wbc][roc] += 100. / max_triggers
            print '  {:03d}\t{}'.format(wbc, '\t'.join(['{:5.1f}%'.format(v) for v in yields[wbc]]))

            # stopping criterion
            if wbc > min_wbc + 3 and any(yld > 10 for yld in yields[wbc - 2]):
                break
        for roc in xrange(self.get_n_rocs()):
            self.api.setDAC('wbc', max(yields, key=lambda x: yields[x][roc]), roc)
        self.api.daqStop()

        # trigger_phase
        print '\nTRIGGER PHASE:'
        for i, trigger_phase in enumerate(trigger_phases):
            if trigger_phase:
                percentage = trigger_phase * 100 / sum(trigger_phases)
                print '{i}\t{d} {v:2.1f}%'.format(i=i, d=int(round(percentage)) * '|', v=percentage)

        self.plot_wbc(yields, plot)

    def plot_wbc(self, yields, show=True):
        if show:
            try:
                x_min = next(key for key, value in yields.items() if any(v > 0 for v in value)) - 1
                x_max = next(key for key, value in OrderedDict(reversed(yields.items())).iteritems() if any(v > 0 for v in value)) + 2
            except StopIteration:
                print('all zero...')
                return
            x = arange(x_min, x_max)
            g = [self.Draw.graph(x, [yields[wbc][roc] for wbc in x], xtit='wbc', ytit='yield [%]', show=False) for roc in range(self.NRocs)]
            self.Draw.multigraph(g, 'WBC Scan', ['ROC {}'.format(i) for i in range(self.NRocs)], 'lp')

    def hitmap(self, t=1, wbc=93, n=None, random_trigger=False):
        pix_data = [pix for event in self.take_data(wbc, t, n, random_trigger) for pix in event.pixels]
        self.Draw.distribution([px.value for px in pix_data], make_bins(-256, 256), xtit='Pulse Height [adc]')
        self.plot_map(pix_data, 'Hit Map', count=True, stats=set_statbox(entries=True))

    def hitmap_random(self, t, n=10000):
        return self.hitmap(t, random_trigger=True, n=n)

    def hitmap_trigger(self, t, wbc=123):
        return self.hitmap(t, wbc, random_trigger=False)

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
            elif data1[0] == 'pix':
                i2c = int(data1[1])
                col = int(data1[2])
                row = int(data1[3])
                self.api.maskPixel(col, row, True, i2c)
                print 'Mask pix', col, row, 'i2c', i2c
            elif data1[0] == 'col':
                i2c = int(data1[1])
                col = int(data1[2])
                print 'Mask col', col, 'i2c', i2c
                for row in xrange(80):
                    self.api.maskPixel(col, row, True, i2c)
            elif data1[0] == 'row':
                i2c = int(data1[1])
                row = int(data1[2])
                print 'Mask row', row, 'i2c', i2c
                for col in xrange(52):
                    self.api.maskPixel(col, row, True, i2c)

    def save_time(self, t=2, n=10000):
        self.api.HVon()
        w = HDF5Writer('main')
        info('taking data ...')
        t_start = time()
        w.PBar.start(t * 60 * 10)
        self.enable_single_pixel(14, 14, prnt=False)
        self.daq_start()
        while time() - t_start < t * 60:
            self.daq_trigger(n)
            for event in self.api.daqGetEventBuffer():
                w.add_event(event)
            sleep(.5)
            w.PBar.update(int((time() - t_start) * 10))
        self.daq_stop()
        w.convert()
        self.api.HVoff()

    def save_random(self, n=10000, n_pixel=10):
        self.api.HVon()
        w = HDF5Writer('main')
        info('taking data ...')
        w.PBar.start(n_pixel)
        for i in xrange(n_pixel):
            self.enable_single_pixel(randint(0, 52), randint(0, 80), prnt=False)
            self.daq_start()
            self.daq_trigger(n)
            for event in self.api.daqGetEventBuffer():
                w.add_event(event)
            self.daq_stop()
            w.PBar.update(i)
        w.convert()
        self.api.HVoff()

    def save_hdf5(self, t=1, n=None, random=False):
        w = HDF5Writer('main')
        self.enable_all()
        w.add_data(self.take_data(w.WBC, t, n, random))
        w.convert()

    def save_data(self, n=240000):
        global BREAK
        t = TreeWriterLjubljana()
        info('START DATA ACQUISITION FOR RUN {}'.format(t.RunNumber))
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
                if i == n:
                    call('ssh -tY f9pc DISPLAY=:0 /home/f9pc001/miniconda2/bin/python /home/f9pc001/Downloads/run/say.py'.split() + ['"finished run {}"'.format(t.RunNumber)])
            except RuntimeError:
                pass
            if BREAK:
                break
        self.daq_stop()
        BREAK = False

    def setup_analogue(self, target_ia=24):
        info('adjusting vana to target {} mA ...'.format(target_ia))
        self.find_vana(target_ia)
        info('adjusting clock delays ...')
        self.find_clk_delay()
        info('adjusting sampling delays ...')
        self.find_tb_delays()
        # todo add caldel scan
        info('adjusting decoding offsets ...')
        self.find_offsets()

    def find_tb_delays(self):
        """findAnalogueTBDelays: configures tindelay and toutdelay"""
        self.set_tb_delays({'tindelay': 0, 'toutdelay': 20})
        self.enable_all()
        self.api.maskAllPixels(1)
        data = mean(self.get_raw_data(), axis=0)
        tin = argmin(data)
        tout = 20 - (data.size - tin - 3 * self.NRocs)
        self.set_tb_delays({'tindelay': tin, 'toutdelay': tout})
        self.save_tb_delays()
        self.api.maskAllPixels(0)

    def find_offsets(self, n_trig=1000):
        b_off, l1_off, alpha = [], [], []
        self.enable_single_pixel(15, 59, prnt=False)
        self.get_raw_event()  # first reading may be bad
        data = self.get_raw_event(n_trig=n_trig)
        if data.shape[1] != 9 * self.get_n_rocs():
            warning('corrupt data!')
            return
        for roc in xrange(self.get_n_rocs()):
            l1_off.append(mean(data[:, (3 + 9 * roc):(8 + 9 * roc)]))
            b_off.append(mean(data[:, 1 + 9 * roc]))
        self.enable_single_pixel(21, 5, prnt=False)
        d_alpha = self.get_raw_event(n_trig=n_trig)
        for roc in range(self.get_n_rocs()):
            high, low0, low1 = mean(d_alpha[:, 5]), mean(data[:, 6]), mean(d_alpha[:, 6])
            alpha.append(1. / (high / (low1 - low0) + 1.))
        self.api.setDecodingL1Offsets(l1_off)
        self.api.setBlackOffsets(b_off)
        self.api.setDecodingAlphas(alpha)
        self.save_config('l1Offset', '[{}]'.format(','.join('{:1.1f}'.format(v) for v in l1_off)))
        self.save_config('blackOffset',  '[{}]'.format(','.join('{:1.1f}'.format(v) for v in b_off)))
        self.save_config('alphas',  '[{}]'.format(','.join('{:1.2f}'.format(v) for v in alpha)))
        self.enable_all()

    def find_clk_delay(self, xmin=0, xmax=20):
        # todo: extend for several rocs
        tin, tout = [self.TBDelays[i] for i in ['tindelay', 'toutdelay']]
        self.set_tb_delays({'tindelay': 0, 'toutdelay': 20})
        x = arange(xmin, xmax)
        self.PBar.start(x.size, counter=True)
        ub, b, ld = array([self.get_mean_header(clk) for clk in x]).T
        x, ub, b, ld = x[ub != 0], ub[ub != 0], b[ub != 0], ld[ub != 0]
        # clk_b = x[argmin(abs(b))]
        # self.Draw.multigraph([self.Draw.graph(x, y, show=False) for y in [ub, b, ld]], 'clk scan', ['ub', 'b', 'ld'])
        self.set_clock(int(round(mean([x[argmax(ld)], mean(x[where(ub < min(ub) + 5)])]))))
        self.set_tb_delays({'tindelay': tin, 'toutdelay': tout})
        self.save_tb_delays()

    def find_vana(self, target=24, xmin=60, xmax=180):
        values = []
        for roc in range(self.NRocs):
            self.set_dac('vana', 0)  # set vana of all ROCs to 0
            values.append(self._find_vana(target, vana=int(mean([xmin, xmax])), step=(xmax - xmin) // 2, roc=roc))
        for i in range(self.NRocs):
            self.set_dac('vana', values[i], i)

    def _find_vana(self, target=24, vana=120, step=60, roc=0):
        self.set_dac('vana', vana, roc)
        c = self.get_ia(prnt=False)
        if abs(c - target) < .2:
            info('analogue current: {:.1f} mA'.format(c))
            self.save_dac_parameters(roc)
            return vana
        return self._find_vana(target, vana + step // 2 * int(sign(target - c)), max(step // 2, 1))

    def draw_address_levels(self, n_trigger=1000, **kwargs):
        x = self.get_address_levels(n_trigger).flatten()
        return self.Draw.distribution(x, make_bins(-512, 512), x_tit='Level [adc]', stats=set_statbox(entries=True), **kwargs)

    def s_curve(self, col=14, row=14, ntrig=1000):
        """ checkADCTimeConstant [vcal=200] [ntrig=10]: sends an amount of triggers for a fixed vcal in high/low region and prints adc values"""
        self.enable_single_pixel(row, col)
        efficiencies = [0 if not px else px[0].value / ntrig for px in self.api.getEfficiencyVsDAC('vcal', 1, 0, 255, nTriggers=ntrig)]
        g = self.Draw.make_tgrapherrors('gsc', 'S-Curve for Pixel {} {}'.format(col, row), x=arange(256), y=efficiencies)
        format_histo(g, x_tit='VCAL', y_tit='Efficiency [%]', y_off=1.3)
        self.Draw.histo(g, draw_opt='ap', lm=.12)


if __name__ == '__main__':
    # command line argument parsing

    parser = ArgumentParser(prog=prog_name, description="A Simple Command Line Interface to the pxar API.")
    parser.add_argument('dir', metavar="DIR", help="The data directory with all required config files. [default = . ]", nargs='?', default='.')
    parser.add_argument('--run', '-r', metavar="FILE", help="Load a cmdline script to be executed before entering the prompt.", default='')
    parser.add_argument('--verbosity', '-v', metavar="LEVEL", default="INFO", help="The output verbosity set in the pxar API.")
    parser.add_argument('--trim', '-T', nargs='?', default=None, help="The output verbosity set in the pxar API. [default = None]")
    parser.add_argument('-wbc', action='store_true')
    args = parser.parse_args(argv)

    print_banner('# STARTING ipython pXar Command Line Interface')

    # start command line
    d = '.' if not args.dir else args.dir
    z = CLIX(d, args.verbosity, args.trim)
    if args.wbc:
        z.wbc_scan()
        while not raw_input('Press enter to restart, anything else to stop: '):
            z.wbc_scan()

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
    ev = z.get_event
    raw = z.daq_get_raw_event
    dt = z.daq_trigger
    # sd = z.api.setDAC
