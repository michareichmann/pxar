# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on June 19th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from os.path import isfile, exists
from os import makedirs, _exit
from ConfigParser import ConfigParser
from datetime import datetime
from ROOT import TFile, gROOT
from numpy import average, sqrt, array
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar


type_dict = {'int32': 'I',
             'uint16': 's',
             'float64': 'D',
             'int64': 'L'}


GREEN = '\033[92m'
ENDC = '\033[0m'
YELLOW = '\033[93m'
RED = '\033[91m'

ON = True
OFF = False


def get_t_str():
    return datetime.now().strftime('%H:%M:%S')


def info(msg, overlay=False, prnt=True):
    if prnt:
        print '{ov}{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head='{}INFO:{}'.format(GREEN, ENDC), ov='\033[1A\r' if overlay else '')


def warning(msg):
    print '{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head='{}WARNING:{}'.format(YELLOW, ENDC))


def critical(msg):
    print '{head} {t} --> {msg}\n'.format(t=get_t_str(), msg=msg, head='{}CRITICAL:{}'.format(RED, ENDC))
    _exit(1)


def file_exists(filename):
    return isfile(filename)


def round_down_to(num, val):
    return int(num) / val * val


def ensure_dir(path):
    if not exists(path):
        info('Creating directory: {d}'.format(d=path))
        makedirs(path)


def is_num(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def choose(v, default, decider='None', *args, **kwargs):
    use_default = decider is None if decider != 'None' else v is None
    if callable(default) and use_default:
        default = default(*args, **kwargs)
    return default if use_default else v


def make_list(value, dtype=None):
    v = array([choose(value, [])]).flatten()
    return v.tolist() if dtype == list else v.astype(dtype) if dtype is not None else v


def print_banner(msg, symbol='=', new_lines=True):
    print '{n}{delim}\n{msg}\n{delim}{n}'.format(delim=(len(str(msg)) + 10) * symbol, msg=msg, n='\n' if new_lines else '')


def do_nothing():
    pass


def load_config(name, ext='ini'):
    parser = ConfigParser()
    parser.read('{}.{}'.format(name, ext))
    return parser


def has_root():
    try:
        import ROOT
        return True
    except ImportError:
        return False


def read_root_file(filename):
    if file_exists(filename):
        return TFile(filename)
    critical('The file: "{}" does not exist...'.format(filename))


def do(fs, pars, exe=-1):
    fs, pars = ([fs], [pars]) if type(fs) is not list else (fs, pars)
    exe = pars if exe == -1 else [exe]
    for f, p, e in zip(fs, pars, exe):
        f(p) if e is not None else do_nothing()


def calculate_column(c0, c1, r0):
    return 2 * (6 * c1 + c0) + (r0 & 1)


def calculate_row(r0, r1, r2):
    return 80 - (36 * r2 + 6 * r1 + r0) / 2


def calculate_col_row(c1, c0, r2, r1, r0):
    return calculate_column(c0, c1, r0), calculate_row(r0, r1, r2)


def bit_shift(value, shift):
    return (value >> shift) & 0b0111


def mean_sigma(values, weights=None):
    """ Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape. """
    weights = [1] * len(values) if weights is None else weights
    if all(weight == 0 for weight in weights):
        return [0, 0]
    avrg = average(values, weights=weights)
    variance = average((values - avrg) ** 2, weights=weights)  # Fast and numerically precise
    return avrg, sqrt(variance)


def set_root_warnings(status):
    gROOT.ProcessLine('gErrorIgnoreLevel = {e};'.format(e='0' if status else 'kError'))


def set_root_output(status=True):
    gROOT.SetBatch(not status)
    set_root_warnings(status)


def remove_letters(string):
    return filter(lambda x: x.isdigit(), string)


def remove_digits(string):
    return filter(lambda x: not x.isdigit(), string)


class PBar:
    def __init__(self):
        self.PBar = None
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.LastUpdate = 0

    def start(self, n):
        self.PBar = ProgressBar(widgets=self.Widgets, maxval=n).start()

    def update(self, i):
        self.finish() if i >= self.PBar.maxval - 1 else self.PBar.update(i + 1)

    def finish(self):
        self.PBar.finish()
