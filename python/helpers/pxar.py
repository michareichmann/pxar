#!/usr/bin/env python
""" Helper classes and functions useful when interfacing the pxar API with Python. """

from lib.PyPxarCore import PixelConfig, PyPxarCore, Statistics, PyRegisterDictionary
from os.path import join, isfile, isdir, basename, realpath
from datetime import datetime
from configparser import ConfigParser
from json import loads
from helpers.utils import info, critical, choose, warning
from numpy import full, array, arange, genfromtxt


NCols = 52
NRows = 80
Frequency = 40e6


# -----------------------------------------
# region CONFIG FILES
class Config(dict):
    """Parent class for the pxar configurations"""
    def __init__(self, filename):
        super().__init__()
        self.FileName = filename
        self.Lines = self.read()

    def read(self, start_at=0, dtype=str):
        if not isfile(self.FileName):
            critical(f'{self.FileName} does not exist!')
        with open(self.FileName) as f:
            lines = f.readlines()
            for words in [line.strip(' \n').split() for line in lines if not line.startswith('--') and not line.startswith('#')]:
                if len(words) > 1:
                    self[words[start_at].lower()] = dtype(' '.join(words[start_at + 1:]))
            return lines

    def show(self, prnt_file=False, prnt=True):
        i = len(max(self, key=len))
        if prnt:
            if prnt_file:
                print(basename(self.FileName))
            for key, value in self.items():
                print(f'{key:{f" <{i}"}} = {value}')

    def get(self, option, default=None):
        return super(Config, self).get(option.lower(), default)

    def get_int(self, option, default=None):
        return int(self.get(option, default))

    def get_b(self, option, default=None):
        return self.get(option, default).encode()

    def add_line(self, key, value):
        self.Lines.append(f'{key} {value}\n')

    def change_line(self, key, value):
        i = next(i for i, line in enumerate(self.Lines) if ' ' in line and key.lower() == line.split()[0].lower())
        words = self.Lines[i].split(' ')
        words[1] = str(value)
        self.Lines[i] = f'{" ".join(words)}\n'

    def set(self, key, value):
        if key is not None:
            if key.lower() in self and self[key.lower()] != value:
                info(f'{basename(self.FileName)}: overwriting "{key}" from {self[key.lower()]} -> {value}')
                self.change_line(key, value)
            elif key.lower() not in self:
                info(f'{basename(self.FileName)}: adding "{key}" key with value {value}')
                self.add_line(key, value)
            self[key.lower()] = value

    @property
    def b(self):
        return {key.encode(): value for key, value in self.items()}


class PxarConfig(Config):
    """ class that loads the old-style config files of psi46expert """
    def __init__(self, filename):
        super().__init__(filename)

    def get_roc_vector(self, opt, n, default=None):
        value = self.get(opt, default)
        return loads(value) if '[' in str(value) else full(n, float(value))

    def save(self, key=None, value=None):
        self.set(key, value)
        with open(join(self.FileName), 'r+') as f:
            f.writelines(self.Lines)
            f.truncate()


class PxarParameters(Config):
    """ class that loads the old-style parameters files of psi46expert """
    AllDACs = PyRegisterDictionary().getAllROCNames()

    def read(self, start_at=1, dtype=int):
        super().read(start_at, dtype)

    def set(self, dac, value, prnt=True, name='ROC'):
        """sets the value of the DAC [dac] to [value]. :returns: old value"""
        if dac.encode() not in self.AllDACs:
            return warning(f'The {name} DAC "{dac}" does not exist!')
        info((f'setting "{dac}" to {value}' if int(value) != self[dac.lower()] else f'{dac} already has a value of {value}'), prnt=prnt)
        old = self[dac]
        self[dac] = int(value)
        return old

    def save(self):
        old_data = genfromtxt(self.FileName, dtype=['i', 'U30', 'i'])
        for i, (_, key, value) in enumerate(old_data):
            if key in self and self[key] != value:
                info(f'{basename(self.FileName)}: overwriting "{key}" from {value} -> {self[key]}')
                old_data[i][-1] = self[key]
        s = [max(len(str(i)) for i in old_data[a]) for a in ['f0', 'f1', 'f2']]
        with open(self.FileName, 'r+') as f:
            f.writelines(['{0:>{3}}  {1:>{4}}  {2:>{5}}\n'.format(*list(tup) + s) for tup in old_data])
            f.truncate()


class TBParameters(PxarParameters):
    """ class that loads the old-style parameters files of psi46expert """
    AllDACs = PyRegisterDictionary().getAllDTBNames()

    def set(self, dac, value, prnt=True, **kwargs):
        return super(TBParameters, self).set(dac, value, prnt, 'testboard parameter')


def create_pixel(col, row, trim=15, roc=0, mask=True):
    p = PixelConfig(col, row, trim)
    p.roc = roc
    p.mask = p in mask if isinstance(mask, PxarMaskFile) else mask
    return p


class PxarMaskFile(list):
    """ class that loads the mask files of pxarGUI """
    def __init__(self, filename):
        super().__init__()
        self.FileName = filename
        self.read()

    def read(self):
        with open(self.FileName) as f:
            lines = [line.strip(' \n').split() for line in f.readlines() if not line.startswith("--") and not line.startswith("#")]
            for words in lines:
                id_, roc, values = words[0], int(words[1]), [int(word) for word in words[2:]]
                words[1:] = [int(word) for word in words[1:]]
                if len(words) == 4:  # single pixel to be masked:
                    self.append(create_pixel(*values, roc=roc))
                elif words[0] == 'col':  # full column to be masked:
                    self.extend([create_pixel(values[0], row, roc=roc) for row in range(NRows)])
                elif words[0] == 'row':  # full row to be masked:
                    self.extend([create_pixel(col, values[0], roc=roc) for col in range(NCols)])
                elif len(words) == 2:  # full roc to be masked
                    self.extend([create_pixel(col, row, roc=roc) for col in range(NCols) for row in range(NRows)])


class PxarTrimConfig(list):
    """ class that loads the old-style trim parameters files of psi46expert """
    def __init__(self, filename, roc, mask):
        super().__init__()
        for vals in genfromtxt(filename, usecols=[2, 3, 0], dtype='i2'):
            self.append(create_pixel(*vals, roc=roc, mask=mask))
# endregion CONFIG FILES
# -----------------------------------------


class PxarStatistics(dict):

    def __init__(self, channels=1):
        super().__init__()
        self.NChannels = max(channels, 1)
        self['General Information'] = {key: 0 for key in ['info words read', 'empty events', 'valid events', 'valid pixels']}
        self['Errors Event'] = {key: 0 for key in ['start', 'stop', 'overflow', 'invalid words', 'invalid XOR', 'frame', 'idledata', 'nodata', 'PKAM']}
        self['Errors TBM'] = {key: 0 for key in ['header', 'trailer', 'eventid mismatch']}
        self['Errors ROC'] = {key: 0 for key in ['missing', 'readback']}
        self['Errors Pixel'] = {key: 0 for key in ['incomplete', 'address', 'pulseheight', 'buffer corrupt']}

    def __str__(self):
        return '\n'.join(f'{head}\n' + '\n'.join([f'    {f"{key}:": <18} {value}' for key, value in dic.items()]) for head, dic in self.items())

    def save(self, hv=None, cur=None):
        hv, cur = choose(f'-{hv}', '', hv), choose(f'-{cur}', '', cur)
        with open(f'stats{hv}{cur}_{datetime.now():%m-%d_%H_%M_%S}.ini', 'w') as f:
            p = ConfigParser()
            p.read_dict(self)
            p.write(f)
            info(f'saved stats in: {f.name}')

    def add(self, stats: Statistics):
        for head, dic in self.items():
            for key in dic:
                self[head][key] += getattr(stats, '_'.join([head.lower()] if 'Errors' in head else [] + [key.lower()]).replace(' ', '_'))

    def clear(self):
        for dic in self.values():
            for key in dic:
                dic[key] = 0

    @property
    def valid_pixels(self):
        return self['General Information']['valid pixels']

    @property
    def valid_events(self):
        return self['General Information']['valid events']

    @property
    def total_events(self):
        return self.valid_pixels + self['General Information']['empty events']

    @property
    def event_rate(self):
        return self.valid_events * Frequency / (self.total_events / float(self.NChannels))

    @property
    def hit_rate(self):
        return self.valid_pixels * Frequency / (self.total_events / float(self.NChannels))


class PxarStartUp:
    """ Initialises the pxar API """
    def __init__(self, directory='.', verbosity='INFO', trim=''):

        # INFO
        self.Dir = realpath(directory) if isdir(directory) else critical('Error: no or invalid configuration directory specified!')
        self.Verbosity = verbosity
        self.Trim = trim
        self.Config = PxarConfig(join(directory, 'configParameters.dat'))
        self.TestBoardName = self.Config.get('testboardName')
        self.ROCType = self.Config.get('rocType')
        self.NROCs, self.I2Cs = self.find_i2cs()
        self.Mask = PxarMaskFile(join(directory, self.Config.get('maskfile')))
        self.IsAnalogue = self.ROCType == 'psi46v2'

        # DACs
        self.TBParameters = self.init_tb_parameters()
        self.TBMDacs = self.init_tbm()
        self.ROCDACs = self.init_dacs()
        self.TrimDACs = self.init_trim_dacs()

        # SETTINGS
        self.PowerSettings = self.init_power()
        self.PGSetup = self.init_pattern_generator()
        self.HubIDs = [int(i) for i in self.Config.get('hubId', 31).split(',')]

        self.API = self.init_api()
        self.set_probes()
        self.set_decoding_offsets()
        info('pxar API is now started and configured.')

    def restart_api(self):
        self.API = self.init_api()
        self.set_probes()
        self.set_decoding_offsets()

    # -----------------------------------------
    # region INIT
    def find_i2cs(self):
        words = self.Config.get('nRocs').split('i2c:')
        n, i2cs = int(words[0]), arange(int(words[0]), dtype='u2') if len(words) == 1 else array(words[1].split(','), 'u2')
        info(f'Number of ROCs:  {n}')
        info(f'Configured I2Cs: {i2cs}')
        return n, i2cs

    def init_tb_parameters(self):
        pars = TBParameters(join(self.Dir, self.Config.get('tbParameters')))
        info(f'tindelay/toutdelay: {pars.get("tindelay")}/{pars.get("toutdelay")}') if self.IsAnalogue else f'clk/phase: {pars.get("clk")}/{pars.get("deser160phase")}'
        return pars

    def init_power(self):
        """ Initialise the power settings of the testboard. """
        power_settings = {key: float(self.Config.get(key, default)) for key, default in [(b'va', 1.9), (b'vd', 2.6), (b'ia', 1.190), (b'id', 1.10)]}
        if any(value > 100 for value in power_settings.values()):
            info('set power settings from [mV] to [V]')
            power_settings = {key: int(value) / 1000. for key, value in power_settings.items()}
        return power_settings

    def init_tbm(self):
        """ Initialise the DACs of the TBMs (TokenBitManager). """
        dacs = [PxarConfig(join(self.Dir, f'{self.Config.get("tbmParameters")}_C{tbm}{i}.dat')) for tbm in range(self.Config.get_int('nTbms')) for i in ['a', 'b']]
        info(f'Found DAC config for {len(dacs)} TBM cores', prnt=len(dacs))
        for pars in dacs:
            pars.show(prnt_file=True, prnt=self.Verbosity != 'INFO')
        return dacs

    def init_dacs(self):
        return [PxarParameters(join(self.Dir, f'{self.Config.get("dacParameters")}{self.Trim}_C{i2c}.dat')) for i2c in self.I2Cs]

    def init_trim_dacs(self):
        dacs = [PxarTrimConfig(join(self.Dir, f'{self.Config.get("trimParameters")}{self.Trim}_C{i2c}.dat'), i2c, self.Mask) for i2c in self.I2Cs]
        s = array([len(i) for i in dacs])
        info(f'There are {s[0]} pixels for all {self.NROCs} ROCs' if all(s == s[0]) else f'ROC pixels: {s}')
        return dacs

    def init_pattern_generator(self):
        pgcal = self.ROCDACs[0].get('wbc') + (6 if 'dig' in self.ROCType else 5)
        info(f'PGCAL: {pgcal}')
        if len(self.TBMDacs) == 0:  # Pattern Generator for single ROC operation:
            return (b'PG_RESR', 25), (b'PG_CAL', pgcal),  (b'PG_TRG', 16), (b'PG_TOK', 0)
        return (b'PG_RESR', 15), (b'PG_CAL', pgcal), (b'PG_TRG', 0)

    def init_api(self):
        try:
            api = PyPxarCore(usbId=self.TestBoardName.encode(), logLevel=self.Verbosity.encode())
            info(f'Init API version: {api.getVersion()}')
            if not api.initTestboard(self.TBParameters.b, self.PowerSettings, self.PGSetup):
                info('Please check if a new FW version is available')
                critical('could not init DTB -- possible firmware mismatch.')
            api.initDUT(self.HubIDs, self.Config.get_b('tbmType', 'tbm08'), byte_dic(self.TBMDacs), self.ROCType.encode(), byte_dic(self.ROCDACs), self.TrimDACs, self.I2Cs)
            api.testAllPixels(True)
            return api
        except RuntimeError as e:
            critical(e)
    # endregion INIT
    # -----------------------------------------

    def set_decoding_offsets(self):
        if not any(word in self.ROCType for word in ['dig', 'proc']):
            b, l1, ub = [self.Config.get_roc_vector(opt, self.NROCs, 0) for opt in ['blackOffset', 'l1Offset', 'alphas']]
            info('set analogue decoding offset set to: {b}')
            info('set analogue level1 offset set to: {l1}')
            info('set analogue alphas set to: {ub}')
            self.API.setBlackOffsets(b)
            self.API.setDecodingL1Offsets(l1)
            self.API.setDecodingAlphas(ub)

    def set_probes(self):
        self.API.SignalProbe(b'a1', self.Config.get_b('probeA1', 'sdata1'))
        self.API.SignalProbe(b'a2', self.Config.get_b('probeA2', 'sdata2'))
        self.API.SignalProbe(b'd1', self.Config.get_b('probeD1', 'clk'))
        self.API.SignalProbe(b'd2', self.Config.get_b('probeD2', 'ctr'))

    def save_dac_parameters(self, roc_id=0):
        return [self.save_dac_parameters(i) for i in range(self.NROCs)] if roc_id is None else self.ROCDACs[roc_id].save()


def byte_dic(lst):
    return [dic.b for dic in lst]
