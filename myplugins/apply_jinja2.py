"""
ApplyJinja2 Plugin for Pelican
-------

This plugin provides a filter to apply jinja2.
"""

from pelican import signals
import jinja2

import logging
logger = logging.getLogger(__name__)


variables = dict()
def configure(generators):
    for generator in generators:
        variables.update(generator.context)

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
