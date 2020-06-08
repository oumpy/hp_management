"""
ApplyJinja2 Plugin for Pelican
-------

This plugin provides a filter to apply jinja2.
"""

from pelican import signals
from pelican.generators import ArticlesGenerator, PagesGenerator
import jinja2

import logging
logger = logging.getLogger(__name__)


variables = dict()
metadata_field = 'jinja2'
def configure(generators):
    for generator in generators:
        variables.update(generator.context)
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            for article in generator.articles:
                if metadata_field in article.metadata and bool(article.metadata[metadata_field]):
                    env = jinja2.Environment(loader=jinja2.DictLoader({'content': article.content, 'title': article.title}))
                    env.filters.update(variables['JINJA_FILTERS'])
                    article._content = env.get_template('content').render(**variables)
                    article.title = env.get_template('title').render(**variables)
        elif isinstance(generator, PagesGenerator):
            for page in generator.pages:
                if metadata_field in page.metadata and bool(page.metadata[metadata_field]):
                    env = jinja2.Environment(loader=jinja2.DictLoader({'content': page.content, 'title':page.title}))
                    env.filters.update(variables['JINJA_FILTERS'])
                    page._content = env.get_template('content').render(**variables)
                    page.title = env.get_template('title').render(**variables)

def filter_apply_jinja2(content):
    env = jinja2.Environment(loader=jinja2.DictLoader({'content': content}))
    env.filters.update(variables['JINJA_FILTERS'])
    template = env.get_template('content')
    result = template.render(**variables)
    return result

def initialize(pelicanobj):
    settings = pelicanobj.settings
    if not 'JINJA_FILTERS' in settings:
        settings['JINJA_FILTERS'] = dict()
    settings['JINJA_FILTERS']['apply_jinja2'] = filter_apply_jinja2
    settings['FILTER_APPLY_JINJA2'] = True

def register():
    signals.initialized.connect(initialize)
    signals.all_generators_finalized.connect(configure)
