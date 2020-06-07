"""
SASS Plugin for Pelican
-------

This plugin converts .scss files in static/sass/ to .css files in static/css/.
"""

from pelican import signals
import sass
import os, shutil
import requests
import zipfile

import logging
logger = logging.getLogger(__name__)

def run_plugin(pelicanobj):
    THEME = pelicanobj.settings['THEME']
    BOOTSTRAP_VERSION = pelicanobj.settings['BOOTSTRAP_VERSION']
    SOURCE_URL = 'https://github.com/twbs/bootstrap/archive/v{}.zip'.format(BOOTSTRAP_VERSION)
    bootstrap_dir = '{}/bootstrap'.format(THEME)
    if not os.path.isdir(bootstrap_dir):
        response = requests.get(SOURCE_URL)
        bootstrapzip = '{}/v{}.zip'.format(THEME, BOOTSTRAP_VERSION)
        with open(bootstrapzip, 'wb') as f:
            f.write(response.content)
        with zipfile.ZipFile(bootstrapzip) as zf:
            zf.extractall(THEME)
        os.remove(bootstrapzip)
        shutil.move('{}/bootstrap-{}'.format(THEME,BOOTSTRAP_VERSION),bootstrap_dir)

    sass.compile(
        dirname=(
            '{}/static/sass'.format(THEME),
            '{}/static/css'.format(THEME),
        ),
        output_style='compressed',
    )

def register():
    signals.initialized.connect(run_plugin)
