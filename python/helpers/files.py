from helpers.utils import make_byte_string
from time import sleep
from glob import glob
from os import rename, remove
from os.path import getsize


def remove_file(name):
    print(f'removing file: {name} ({make_byte_string(getsize(name))})')
    remove(name)


def remove_trim_files(i2c=None):
    for name in glob(f'*Parameters*{"" if i2c is None else f"C{i2c}"}.dat'):
        if '_' in name and name.split('_')[-2][-1].isdigit():
            remove_file(name)


def remove_logs():
    for name in glob('pxar*'):
        if '.root' in name or '.log' in name or name.endswith('~'):
            remove_file(name)


def remove_swaps():
    for name in glob('*~'):
        remove_file(name)


def get_old_i2c(i2c=None):
    i2cs = [''.join(filter(str.isdigit, name)) for name in glob('dacParameters_*.dat')]
    if i2c and str(i2c) not in i2cs:
        raise ValueError(f'I²C {i2c} not found!')
    return min(i2cs, key=int) if i2c is None else str(i2c)


def rename_files(new, old):
    remove_trim_files()  # if the i2c changed the files with the trimmed parameters are not required
    for name in glob(f'*Parameters_C{old}.dat'):
        print(f'renaming {name} --> {name.replace(old, new)}')
        rename(name, name.replace(old, new))


def change_config(new, old):
    with open('configParameters.dat', 'r+') as f:
        lines = f.readlines()
        i = next(i for i, line in enumerate(lines) if line.startswith('nRocs'))
        old_i2c_str = lines[i].split(":")[1].strip(' \n')
        new_i2c_str = ','.join(sorted(old_i2c_str.replace(old, new).split(','), key=int))
        print(f'overwriting I²C in configParameters.dat: {old_i2c_str} --> {new_i2c_str}')
        lines[i] = ':'.join([lines[i].split(':')[0], f' {new_i2c_str}\n'])
        f.seek(0)
        f.writelines(lines)


def find_i2c(max_i2c=16):
    from helpers.pxar import PxarStartUp
    for i2c in range(max_i2c):
        rename_files(i2c, get_old_i2c())
        pxar = PxarStartUp('.')
        api = pxar.API
        sleep(2)
        if api.getTBia() * 1000 > 5:
            print('found i2c', i2c)
            return True
        print(i2c, api.getTBia() * 1000)
        del api
