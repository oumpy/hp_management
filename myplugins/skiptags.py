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

def run(generator):
    reglist = [re.compile('<{0}.*?>(.*?)</{0}>'.format(tag)) for tag in generator.settings['SKIPTAGS']]
    if isinstance(generator, ArticlesGenerator):
        attr = 'articles'
    else:
        attr = 'pages'
    for obj in getattr(generator, attr):
        for reg in reglist:
            obj._content = re.sub(reg, '', obj.content)

def register():
    signals.initialized.connect(initialize)
    signals.page_generator_finalized.connect(run)
    signals.article_generator_finalized.connect(run)
