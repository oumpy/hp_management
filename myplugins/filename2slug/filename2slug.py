"""
filename2slug Plugin for Pelican
-------

This plugin sets the filename of the source (not title) as the default slug for articles/pages.
"""

from pelican import signals
import os

import logging
logger = logging.getLogger(__name__)


def filename2slug(obj):
    if 'slug' in obj.metadata.keys():
        return
    
    source_path = obj.source_path
    filename = os.path.basename(source_path)
    filename_body = filename[:filename.rfind('.')]
    obj.slug = filename_body


def register():
    signals.content_object_init.connect(filename2slug)
