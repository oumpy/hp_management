"""
path2obj Plugin for Pelican
-------

This plugin makes two dictionary objects `source2obj` & `url2obj`,
returning the corresponding object to a given path (for articles/pages),
or category/tag/author name.
"""

from pelican import signals
from pelican.generators import ArticlesGenerator, StaticGenerator, PagesGenerator

import logging
logger = logging.getLogger(__name__)


source2obj = dict()
url2obj = dict()

def run_plugin(generators):
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            for article in generator.articles:
                source2obj[article.relative_source_path] = article
                url2obj[article.url] = article
            for cat in generator.categories:
                source2obj['category:' + cat[0].name] = cat[0]
                url2obj[cat[0].url] = cat[0]
            for tag in generator.tags:
                source2obj['tag:' + tag.name] = tag
                url2obj[tag.url] = tag
            for author in generator.authors:
                source2obj['author:' + author[0].name] = author[0]
                url2obj[author[0].url] = author[0]
        elif isinstance(generator, PagesGenerator):
            for page in generator.pages:
                source2obj[page.relative_source_path] = page
                url2obj[page.url] = page
        elif isinstance(generator, StaticGenerator):
            for file_ in generator.staticfiles:
                source2obj[file_.relative_source_path] = file_
                url2obj[file_.url] = file_

def filter_source2obj(x):
    return source2obj[x]

def filter_url2obj(x):
    return url2obj[x]

def initialize(pelicanobj):
    pelicanobj.settings['JINJA_FILTERS']['url2obj'] = filter_url2obj
    pelicanobj.settings['JINJA_FILTERS']['source2obj'] = filter_source2obj

def register():
    signals.initialized.connect(initialize)
    signals.all_generators_finalized.connect(run_plugin)
