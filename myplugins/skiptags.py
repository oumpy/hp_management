"""
SkipTags Plugin for Pelican
-------

This plugin skips specified tags in pages/articles.
"""

from pelican import signals
from pelican.generators import ArticlesGenerator, PagesGenerator
import re

import logging
logger = logging.getLogger(__name__)

def initialize(pelicanobj):
    from pelican.settings import DEFAULT_CONFIG
    DEFAULT_CONFIG.setdefault('SKIPTAGS', ['h1',])
    if pelicanobj:
        pelicanobj.settings.setdefault('SKIPTAGS', ['h1',])

def configure(generators):
    for generator in generators:
        reglist = [re.compile('<{0}.*?>(.*?)</{0}>'.format(tag)) for tag in generator.settings['SKIPTAGS']]
        if isinstance(generator, ArticlesGenerator):
            for article in generator.articles:
                for reg in reglist:
                    article._content = re.sub(reg, '', article.content)
        elif isinstance(generator, PagesGenerator):
            for page in generator.pages:
                for reg in reglist:
                    page._content = re.sub(reg, '', page.content)

def register():
    signals.initialized.connect(initialize)
    signals.all_generators_finalized.connect(configure)
