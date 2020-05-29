"""
RecentPosts Plugin for Pelican
-------

This plugin provides a filter to obtain recent articles.
"""

from pelican import signals
from pelican.generators import ArticlesGenerator, PagesGenerator
import datetime

import logging
logger = logging.getLogger(__name__)

def filter_recentposts(articles, days=30, counts=5):
    if len(articles) <= counts:
        return articles
    else:
        last = articles[0].date.toordinal()
        for i in range(len(articles)):
            d = articles[i].date.toordinal()
            if d < last - days:
                break
        else:
            i += 1
        return articles[:max(counts, i)]

def initialize(pelicanobj):
    settings = pelicanobj.settings
    if not 'JINJA_FILTERS' in settings:
        settings['JINJA_FILTERS'] = dict()
    settings['JINJA_FILTERS']['recentposts'] = filter_recentposts
    settings['FILTER_APPLY_JINJA2'] = True

def register():
    signals.initialized.connect(initialize)
