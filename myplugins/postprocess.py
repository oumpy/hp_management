"""
Postprocess Plugin for Pelican
-------

This plugin runs a postprocess command after finishing complile.
"""

from pelican import signals
import subprocess

import logging
logger = logging.getLogger(__name__)


def run_postprocess(pelicanobj):
    command = pelicanobj.settings['POSTPROCESS_COMMAND']
    output = pelicanobj.settings['OUTPUT_PATH']
    subprocess.run(command.split() + [ output ])

def register():
    signals.finalized.connect(run_postprocess)
