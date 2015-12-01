#!/usr/bin/env python2
"""
Helper classes and functions useful when interfacing the pxar API with Python.
"""
import PyPxarCore
from PyPxarCore import Pixel, PixelConfig, PyPxarCore, PyRegisterDictionary, PyProbeDictionary
from functools import wraps # used in parameter verification decorator ("arity")
import os # for file system cmds
import sys

# "arity": decorator used for parameter parsing/verification on each cmd function call
# Usually, the cmd module only passes a single string ('line') with all parameters;
# this decorator divides and verifies the types of each parameter.
def arity(n, m, cs=[]): # n = min number of args, m = max number of args, cs = types
    def __temp1(f):
        @wraps(f) # makes sure the docstring of the orig. function is passed on
        def __temp2(self, text):
            ps = filter(lambda p: p, text.split(" "))
            if len(ps) < n:
                print "Error: this command needs %d arguments (%d given)" % (n, len(ps))
                return
            if len(ps) > m:
                print "Error: this command takes at most %d arguments (%d given)" % (m, len(ps))
                return
            # now verify the type
            try:
                ps = [ c(p) for c, p in zip(cs, ps) ]
            except ValueError as e:
                print "Error: '" + str(p) + "' does not have " + str(c)
                return
            f(self, *ps)
        return __temp2
    return __temp1

def print_data(fullOutput,data,stepsize=1):
    for idac, dac in enumerate(data):
        s = "DAC " + str(idac*stepsize) + ": "
        if fullOutput:
            for px in dac:
                s += str(px)
        else:
            s += str(len(dac)) + " pixels"
        print s

def get_possible_filename_completions(text):
    head, tail = os.path.split(text.strip())
    if head == "": #no head
        head = "."
    files = os.listdir(head)
    return [ f for f in files if f.startswith(tail) ]

def extract_full_argument(line, endidx):
    newstart = line.rfind(" ", 0, endidx)
    return line[newstart:endidx]

class PxarConfigFile:
    """ class that loads the old-style config files of psi46expert """
    def __init__(self, f):
        self.config = {}
        import shlex
        thisf = open(f)
        try:
            for line in thisf:
                if not line.startswith("--") and not line.startswith("#"):
                    parts = shlex.split(line)
                    if len(parts) == 2:
                        self.config[parts[0].lower()] = parts[1]
                    elif len(parts) == 4:
                        parts = [parts[0],' '.join(parts[1:])]
                        if len(parts) == 2:
                            self.config[parts[0].lower()] = parts[1]

        finally:
            thisf.close()
    def show(self):
        print self.config
    def get(self, opt, default = None):
        if default:
            return self.config.get(opt.lower(),default)
        else:
            return self.config[opt.lower()]

class PxarParametersFile:
    """ class that loads the old-style parameters files of psi46expert """
    def __init__(self, f):
        self.config = {}
        import shlex
        thisf = open(f)
        try:
            for line in thisf:
                if not line.startswith("--") and not line.startswith("#"):
                    parts = shlex.split(line)
                    if len(parts) == 3:
                        # ignore the first part (index/line number)
                        self.config[parts[1].lower()] = parts[2]
                    elif len(parts) == 2:
                        self.config[parts[0].lower()] = parts[1]
        finally:
            thisf.close()
    def show(self):
        print self.config
    def get(self, opt, default = None):
        if default:
            return self.config.get(opt.lower(),default)
        else:
            return self.config[opt.lower()]
    def getAll(self):
        return self.config

class TrimFile:
    """class that loads the trim file"""
    def __init__(self, dir, i2c):
        self.file_name = 'trimParameters'
        self.use_trim = True
        self.rows = 52
        self.columns = 80
        self.config = []
        # self.config = self.create_default()

        self.file_path = '{dir}/{name}_C{i2c}.dat'.format(dir=dir, i2c=i2c, name=self.file_name)
        f = open(self.file_path)
        try:
            for line in f:
                data = line.split()
                row = int(data[3])
                col = int(data[2])
                trim = int(data[0])
                conf = PixelConfig(col, row, trim)
                conf.mask = False
                self.config.append(conf)
        except IndexError:
            self.config = self.create_default()
        finally:
            f.close()

        self.parameters = len(self.config)
        self.is_full_file = True if self.parameters == self.rows * self.columns else False

    def show(self):
        print 'col\t row\t trim'
        for i in range(20):
            print self.config[i].column, self.config[i].row, self.config[i].trim
        return

    def create_default(self):
        conf = []
        for col in range(self.columns):
            for row in range(self.rows):
                p = PixelConfig(col, row, 15)
                p.mask = False
                conf.append(p)
        return conf

    def get_trim(self):
        return self.config

    def add_to_default(self):
        default = self.create_default()
        for j, def_conf in enumerate(default):
            for i, conf in enumerate(self.config):
                if def_conf.row == conf.row and def_conf.column == conf.column:
                    default[j].trim = conf.trim
                    self.config.pop(i)
        self.config = default
        return default

def PxarStartup(directory, verbosity):
    if not directory or not os.path.isdir(directory):
        print "Error: no or invalid configuration directory specified!"
        sys.exit(404)

    config = PxarConfigFile('%sconfigParameters.dat'%(os.path.join(directory,"")))
    tbparameters = PxarParametersFile('%s%s'%(os.path.join(directory,""),config.get("tbParameters")))

    # Power settings:
    power_settings = {
        "va":config.get("va",1.9),
        "vd":config.get("vd",2.6),
        "ia":config.get("ia",1.190),
        "id":config.get("id",1.10)}

    tbmDACs = []
    for tbm in range(int(config.get("nTbms"))):
        for n in range(2):
            tbmparameters = PxarParametersFile('%s%s'%(os.path.join(directory,""),config.get("tbmParameters") + "_C" + str(tbm) + ("a" if n%2 == 0 else "b") + ".dat"))
            tbmDACs.append(tbmparameters.getAll())

    print "Have DAC config for " + str(len(tbmDACs)) + " TBM cores:"
    for idx, tbmDAC in enumerate(tbmDACs):
        for key in tbmDAC:
            print "  TBM " + str(idx/2) + ("a" if idx%2 == 0 else "b") + " dac: " + str(key) + " = " + str(tbmDAC[key])

    # init pixel list
    pixels = list()
    for column in range(0, 52):
        for row in range(0, 80):
            p = PixelConfig(column,row,15)
            p.mask = False
            pixels.append(p)

    rocDacs = []
    rocPixels = list()
    rocI2C = []
    config_nrocs = config.get("nrocs").split()
    nrocs = int(config_nrocs[0])
    i2cs = [i for i in range(nrocs)]
    if len(config_nrocs) > 1:
        if config_nrocs[1].startswith('i2c'):
            i2cs = ' '.join(config_nrocs[2:])
            i2cs = [int(i) for i in i2cs.split(',')]
            print "Number of ROCs:", nrocs, "\b; Configured I2C's:", i2cs
    for roc in xrange(nrocs):
        if len(i2cs)> roc:
            i2c = i2cs[roc]
        else:
            i2c = roc
        dacconfig = PxarParametersFile('%s%s_C%i.dat'%(os.path.join(directory,""),config.get("dacParameters"),i2c))
        trimconfig = TrimFile(directory, i2c)
        if trimconfig.use_trim:
            print 'use trimming from file:', trimconfig.file_path
            trimconfig.show()
            if trimconfig.is_full_file:
                pixels = trimconfig.get_trim()
            else:
                pixels = trimconfig.add_to_default()

        rocI2C.append(i2c)
        rocDacs.append(dacconfig.getAll())
        rocPixels.append(pixels)


    # set pgcal according to wbc
    pgcal = int(rocDacs[0]['wbc']) + 6 if nrocs else 106

    # Pattern Generator for single ROC operation:
    if int(config.get("nTbms")) == 0:
        pg_setup = (
            #("PG_RESR",25),
            ("PG_RESR",25),
            ("PG_CAL",pgcal),
            ("PG_TRG",16),
            #("PG_TRG",16),
            ("PG_TOK",0)
            #("PG_RESR",0),
            )
    else:
        pg_setup = (
            ("PG_RESR",15),
            ("PG_CAL",pgcal),
            ("PG_TRG",0))
       # Start an API instance from the core pxar library
    api = PyPxarCore(usbId=config.get("testboardName"),logLevel=verbosity)
    print api.getVersion()
    if not api.initTestboard(pg_setup = pg_setup,
    power_settings = power_settings,
    sig_delays = tbparameters.getAll()):
        print "WARNING: could not init DTB -- possible firmware mismatch."
        print "Please check if a new FW version is available"
        exit()
    try:
        api.setDecodingOffset(int(config.get("decodingOffset", int(0))))
    except KeyError:
        pass
    print "And we have just initialized " + str(len(pixels)) + " pixel configs to be used for every ROC!"

    api.initDUT(int(config.get("hubId",31)),config.get("tbmType","tbm08"),tbmDACs,config.get("rocType"),rocDacs,rocPixels, rocI2C)

    api.testAllPixels(True)
    print "Now enabled all pixels"

    print "pxar API is now started and configured."
    return api
