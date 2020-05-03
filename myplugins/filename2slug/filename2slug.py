"""
filename2slug Plugin for Pelican
-------

This plugin sets the filename of the source (not title) as the default slug for articles/pages.
"""

from __future__ import unicode_literals
from pelican import signals
import os
from pelican.generators import ArticlesGenerator, PagesGenerator

import logging
logger = logging.getLogger(__name__)


def filename2slug(obj, settings):
    if hasattr(obj.metadata, 'slug'):
        return
    
    source_path = obj.source_path
    filename = os.path.basename(source_path)
    filename_body = filename[:filename.rfind('.')]
    obj.slug = filename_body


def run_plugin(generators):
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            for article in generator.articles:
                filename2slug(article, generator.settings)
        elif isinstance(generator, PagesGenerator):
            for page in generator.pages:
                filename2slug(page, generator.settings)


def register():
    signals.all_generators_finalized.connect(run_plugin)
