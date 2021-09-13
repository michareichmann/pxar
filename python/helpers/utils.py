# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on June 19th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import PyConfig
PyConfig.IgnoreCommandLineOptions = True
from os.path import isfile, exists, dirname, realpath
from os import makedirs, _exit
from pickle import loads
from configparser import ConfigParser, NoSectionError, NoOptionError
from datetime import datetime
from time import time
from numpy import average, sqrt, array, count_nonzero, zeros, full, log2, quantile, cos, sin, arctan2
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar, Widget, SimpleProgress
from uncertainties import ufloat
from uncertainties.core import Variable, AffineScalarFunc
from functools import wraps
from copy import deepcopy


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


def prnt_msg(txt, head, color=None, blank_lines=0, endl=True, prnt=True):
    if prnt:
        print('\n' * blank_lines + f'\r{color}{head}:{ENDC} {get_t_str()} --> {txt}', end='\n' if endl else ' ')


def info(txt, blank_lines=0, endl=True, prnt=True):
    prnt_msg(txt, 'INFO', GREEN, blank_lines, endl, prnt)
    return time()


def add_to_info(t, msg='Done', prnt=True):
    if prnt:
        print(('{m} ({t:2.2f} s)'.format(m=msg, t=time() - t)))


def warning(txt, blank_lines=0, prnt=True):
    prnt_msg(txt, 'WARNING', YELLOW, blank_lines, prnt=prnt)


def critical(txt):
    prnt_msg(txt, 'CRITICAL', RED)
    _exit(2)


def print_elapsed_time(start, what='This'):
    print(('Elapsed time for {w}: {t}'.format(t=get_elapsed_time(start), w=what)))


def get_elapsed_time(start):
    t = datetime.fromtimestamp(time() - start)
    return '{}.{:02.0f}'.format(t.strftime('%M:%S'), t.microsecond / 10000)


def round_down_to(num, val=1):
    return int(num) // val * val


def round_up_to(num, val=1):
    return int(num) // val * val + val


def get_base_dir():
    return dirname(dirname(realpath(__file__)))


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


def is_iter(v):
    try:
        iter(v)
        return True
    except TypeError:
        return False


def is_ufloat(value):
    return type(value) in [Variable, AffineScalarFunc]


def make_ufloat(n, s=0):
    if is_iter(n):
        return array([ufloat(*v) for v in array([n, s]).T])
    return n if is_ufloat(n) else ufloat(n, s)


def uarr2n(arr):
    return array([i.n for i in arr]) if is_ufloat(arr[0]) else arr


def print_banner(msg, symbol='=', new_lines=True):
    print('{n}{delim}\n{msg}\n{delim}{n}'.format(delim=(len(str(msg)) + 10) * symbol, msg=msg, n='\n' if new_lines else ''))


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
    if isfile(filename):
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


def mean_sigma(values, weights=None, err=True):
    """ Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape. """
    if len(values) == 1:
        value = make_ufloat(values[0])
        return (value, ufloat(value.s, 0)) if err else (value.n, value.s)
    weights = full(len(values), 1) if weights is None else weights
    if is_ufloat(values[0]):
        errors = array([v.s for v in values])
        weights = full(errors.size, 1) if all(errors == errors[0]) else [1 / e if e else 0 for e in errors]
        values = array([v.n for v in values], 'd')
    if all(weight == 0 for weight in weights):
        return [0, 0]
    avrg = average(values, weights=weights)
    sigma = sqrt(average((values - avrg) ** 2, weights=weights))  # Fast and numerically precise
    m, s = ufloat(avrg, sigma / (sqrt(len(values)) - 1)), ufloat(sigma, sigma / sqrt(2 * len(values)))
    return (m, s) if err else (m.n, s.n)


def remove_letters(string):
    return [x for x in string if x.isdigit()]


def remove_digits(string):
    return [x for x in string if not x.isdigit()]


def calc_eff(k=0, n=0, values=None):
    values = array(values) if values is not None else None
    if n == 0 and (values is None or not values.size):
        return zeros(3)
    k = float(k if values is None else count_nonzero(values))
    n = float(n if values is None else values.size)
    m = (k + 1) / (n + 2)
    mode = k / n
    s = sqrt(((k + 1) / (n + 2) * (k + 2) / (n + 3) - ((k + 1) ** 2) / ((n + 2) ** 2)))
    return array([mode, max(s + (mode - m), 0), max(s - (mode - m), 0)]) * 100


# ----------------------------------------
# region CLASSES
class PBar(object):
    def __init__(self, start=None, counter=False, t=None):
        self.PBar = None
        self.Widgets = self.init_widgets(counter, t)
        self.Step = 0
        self.N = 0
        self.start(start)

    def __reduce__(self):
        return self.__class__, (None, False, None), (self.Widgets, self.Step, self.N)

    def __setstate__(self, state):
        self.Widgets, self.Step, self.N = state
        if self.N:
            self.PBar = ProgressBar(widgets=self.Widgets, maxval=self.N).start()
            self.update(self.Step) if self.Step > 0 else do_nothing()

    @staticmethod
    def init_widgets(counter, t):
        return ['Progress: ', SimpleProgress('/') if counter else Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed() if t is None else EventSpeed(t)]

    def start(self, n, counter=None, t=None):
        if n is not None:
            self.Step = 0
            self.PBar = ProgressBar(widgets=self.Widgets if t is None and counter is None else self.init_widgets(counter, t), maxval=n).start()
            self.N = n

    def update(self, i=None):
        i = self.Step if i is None else i
        if i >= self.PBar.maxval:
            return
        self.PBar.update(i + 1)
        self.Step += 1
        if i == self.PBar.maxval - 1:
            self.finish()

    def finish(self):
        self.PBar.finish()

    def is_finished(self):
        return self.PBar.currval == self.N


class EventSpeed(Widget):
    """Widget for showing the event speed (useful for slow updates)."""

    def __init__(self, t='s'):
        self.unit = t
        self.factor = {'s': 1, 'min': 60, 'h': 60 * 60}[t]

    def update(self, pbar):
        value = 0
        if pbar.seconds_elapsed > 2e-6 and pbar.currval > 2e-6:
            value = pbar.currval / pbar.seconds_elapsed * self.factor
        return '{:4.1f} E/{}'.format(value, self.unit)


def update_pbar(func):
    @wraps(func)
    def my_func(*args, **kwargs):
        if args[0].PBar is not None and args[0].PBar.PBar is not None and not args[0].PBar.is_finished():
            args[0].PBar.update()
        return func(*args, **kwargs)
    return my_func


class Config(ConfigParser):

    def __init__(self, file_name, **kwargs):
        ConfigParser.__init__(self, **kwargs)
        self.FileName = file_name
        self.read(file_name)

    def get_value(self, section, option, dtype=str, default=None):
        dtype = type(default) if default is not None else dtype
        try:
            if dtype is bool:
                return self.getboolean(section, option)
            v = self.get(section, option)
            return loads(v) if dtype == list or '[' in v and dtype is not str else dtype(v)
        except (NoOptionError, NoSectionError):
            return default

    def get_values(self, section):
        return [j for i, j in self.items(section)]

    def get_list(self, section, option, default=None):
        return self.get_value(section, option, list, choose(default, []))

    def show(self):
        for key, section in list(self.items()):
            print(('[{}]'.format(key)))
            for option in section:
                print(('{} = {}'.format(option, self.get(key, option))))
            print()


def make_byte_string(v):
    n = int(log2(v) // 10) if v else 0
    return f'{v / 2 ** (10 * n):1.1f} {["B", "kB", "MB", "GB"][n]}'


def prep_kw(dic, **default):
    d = deepcopy(dic)
    for kw, value in default.items():
        if kw not in d:
            d[kw] = value
    return d


def freedman_diaconis(x):
    return 2 * (quantile(x, .75) - quantile(x, .25)) / x.size ** (1 / 3)


def get_x(x1, x2, y1, y2, y):
    return (x2 - x1) / (y2 - y1) * (y - y1) + x1


def get_y(x1, x2, y1, y2, x):
    return get_x(y1, y2, x1, x2, x)


def cart2pol(x, y):
    return array([sqrt(x ** 2 + y ** 2), arctan2(y, x)])


def pol2cart(rho, phi):
    return array([rho * cos(phi), rho * sin(phi)])
# endregion CLASSES
# ----------------------------------------
