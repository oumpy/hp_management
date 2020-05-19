"""
URL-Encode Plugin for Pelican
-------

This plugin provides a filter to encode URLs.
"""

from pelican import signals
import urllib.parse

import logging
logger = logging.getLogger(__name__)

def filter_url_encode(s):
    return urllib.parse.quote(str(s))

def initialize(pelicanobj):
    settings = pelicanobj.settings
    if not 'JINJA_FILTERS' in settings:
        settings['JINJA_FILTERS'] = dict()
    settings['JINJA_FILTERS']['url_encode'] = filter_url_encode

def register():
    signals.initialized.connect(initialize)
