"""
CategoryNames Plugin for Pelican
-------

This plugin replaces or adds alternative names for categories.
"""

from pelican import signals

import logging
logger = logging.getLogger(__name__)

def initialize(instance):
    from pelican.settings import DEFAULT_CONFIG
    DEFAULT_CONFIG.setdefault('CATEGORYNAMES_ATTRIBUTES',
                              ('name','name_long'))
    DEFAULT_CONFIG.setdefault('CATEGORYNAMES_ALTERNATIVES',
                              dict())
    instance.settings.setdefault('CATEGORYNAMES_ATTRIBUTES',
                                 ('name','name_long'))
    instance.settings.setdefault('CATEGORYNAMES_ALTERNATIVES',
                                 dict())


def rename_category(category, settings):
    attrs = settings['CATEGORYNAMES_ATTRIBUTES']    
    alters = settings['CATEGORYNAMES_ALTERNATIVES']
    cat = str(category)
    category.slug = cat.lower()
    if cat in alters.keys():
        new_cats = alters[cat]
        if not isinstance(new_cats,(list,tuple)):
            new_cats = (new_cats,)
        for i in range(min(len(attrs),len(new_cats))):
            setattr(category, attrs[i], new_cats[i])


def run_plugin(generator):
    for article in generator.articles:
        rename_category(article.category, generator.settings)


def register():
    signals.initialized.connect(initialize)
    signals.article_generator_pretaxonomy.connect(run_plugin)
