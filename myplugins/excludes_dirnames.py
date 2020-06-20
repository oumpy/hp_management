"""
ExcludesDirnames Plugin for Pelican
-------

This plugin provides setting parameters,
ARTICLE_EXCLUDES_DIRNAMES, PAGE_EXCLUDES_DIRNAMES and STATIC_EXCLUDES_DIRNAMES.
"""

from __future__ import unicode_literals
from pelican import signals
import os, pathlib

import logging
logger = logging.getLogger(__name__)

def defaultsettings(pelicanobj):
    from pelican.settings import DEFAULT_CONFIG
    DEFAULT_CONFIG.setdefault('ARTICLE_EXCLUDES_DIRNAMES', [])
    DEFAULT_CONFIG.setdefault('PAGE_EXCLUDES_DIRNAMES', [])
    DEFAULT_CONFIG.setdefault('STATIC_EXCLUDES_DIRNAMES', [])
    if pelicanobj:
        pelicanobj.settings.setdefault('ARTICLE_EXCLUDES_DIRNAMES', [])
        pelicanobj.settings.setdefault('PAGE_EXCLUDES_DIRNAMES', [])
        pelicanobj.settings.setdefault('STATIC_EXCLUDES_DIRNAMES', [])

def run_plugin(pelicanobj):
    defaultsettings(pelicanobj)
    settings = pelicanobj.settings
    PATH = settings['PATH'] if 'PATH' in settings else './'
    variable_sets = [
        ('ARTICLE_EXCLUDES', 'ARTICLE_PATHS', 'ARTICLE_EXCLUDES_DIRNAMES'), 
        ('PAGE_EXCLUDES', 'PAGE_PATHS', 'PAGE_EXCLUDES_DIRNAMES'),
        ('STATIC_EXCLUDES', 'STATIC_PATHS', 'STATIC_EXCLUDES_DIRNAMES'),
    ]

    for excludes, paths, dirnames in variable_sets:
        settings[excludes] += [
            os.path.relpath(dir, PATH)
            for root in settings[paths]
                for dirname in settings[dirnames]
                    for dir in pathlib.Path(PATH).glob(os.path.join(root, '**/{}'.format(dirname)))
                        if dir.is_dir()
        ]

def register():
    signals.initialized.connect(run_plugin)
