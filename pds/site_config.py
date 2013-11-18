#!/usr/bin/env python
"""
site configuration for tdl (borrowed from larch):

"""
from __future__ import print_function

import os
import sys


def get_homedir():
    "determine home directory"
    def check(method, s):
        try:
            if method(s) not in (None, s):
                return method(s)
        except KeyError:
            pass
        return None

    home_dir = check(os.path.expanduser, '~')
    if home_dir is not None:
        for var in ('$HOME', '$USERPROFILE', '$ALLUSERSPROFILE', '$HOMEPATH'):
            home_dir = check(os.path.expandvars, var)
            if home_dir is not None:
                break

    if home_dir is None:
        home_dir = os.path.abspath('.')
    return home_dir.replace('\\', '/')


def make_dir(dname):
    if os.path.exists(dname):
        return True
    try:
        os.mkdir(dname)
        return True
    except (OSError, TypeError):
        print('Error trying to create %s' % dname)
        print(sys.exc_info()[1])
    return False

pds_dirname = '.pds'
if os.name == 'nt':
    pds_dirname = 'pds'

home_dir     = get_homedir()
user_pds_dir = os.path.abspath(os.path.join(home_dir, pds_dirname))
if not os.path.exists(user_pds_dir):
    make_dir(user_pds_dir)

