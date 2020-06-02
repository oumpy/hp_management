"""
MakeMenu Plugin for Pelican
-------

This plugin provides a filter to generate submenus list.
"""

from __future__ import unicode_literals
from pelican import signals
from collections import defaultdict
import re
import sys, os

import logging
logger = logging.getLogger(__name__)

filter_url2obj = None

def initialize(pelicanobj):
    global filter_url2obj
    settings = pelicanobj.settings
    if 'url2obj' in settings['JINJA_FILTERS'].keys():
        filter_url2obj = settings['JINJA_FILTERS']['url2obj']
    settings['JINJA_FILTERS']['makemenu'] = makemenu


class MenuItem():
    def __init__(self, url, title=None, subsections=[], active_pages=None):
        self.url, self.title, self.active_pages = url, title, active_pages
        if subsections == []:
            self.subsections = []
        else:
            self.subsections = subsections
    AUTO = None

def makemenu(add_on_menu, page_url, depth=1, CARET=False):
# add_on_menu : list or tuple of MenuItem objects
# page_url: give the url of the page where this filter is called
# depath: depth of the menu hierarchy
    ret = []
    page_depth = page_url.count('/')
    if page_depth == 0:
        rooturl = '.'
    else:
        rooturl = '../' * (page_depth-1) + '..' 
    FORWARD, BACK = 0, 1
    for obj in add_on_menu:
        pool = [(obj, 0, FORWARD, {'parent':None})]
        active_flag = defaultdict(bool)
        active_flag[page_url] = True
        while pool:
            node, d, s, params = pool.pop()
            if isinstance(node, str):
                if filter_url2obj:
                    node = filter_url2obj(node)
                else:
                    logger.error('You need path2obj for \'{}\' in submenu.'.format(node))
            if not hasattr(node, 'subsections') or node.subsections is None:
                if filter_url2obj:
                    subsections = filter_url2obj(node.url).subsections
                else:
                    subsections = []
            else:
                subsections = node.subsections
            if node.title is None and filter_url2obj:
                title = filter_url2obj(node.url).title
            else:
                title = node.title

            if s == FORWARD:
                if d >= depth:
                    subsections = []

                params['li_line'] = len(ret)
                ret.append('<li class="nav-item">')
                caret = ''
                if d == 0:
                    if subsections and CARET:
                        caret = ' <span class="caret"></span>'
                    a_format = '<a class="nav-link" href="{}">{}{}</a>'
                else:
                    a_format = '<a class="nav-link" href="{}">{}{}</a>'
                ret.append(a_format.format('{}/{}'.format(rooturl, node.url), title, caret))

                if len(subsections) > 0:
                    ret.append('<ul class="navbar-nav mr-auto">')
                    params['/ul'] = True
                    pool.append((node, d, BACK, params))
                    for c in subsections[::-1]:
                        pool.append((c, d+1, FORWARD, {'parent':node.url}))
                else:
                    params['/ul'] = False
                    pool.append((node, d, BACK, params))
            else: # s==BACK
                if params['/ul']:
                    ret.append('</ul></li>')
                else:
                    ret.append('</li>')
                active_flag[node.url] |= bool(
                    hasattr(node,'active_pages') and node.active_pages and re.match(node.active_pages, page_url)
                )
                if active_flag[node.url]:
                    ret[params['li_line']] = ret[params['li_line']][:-2] + ' active">'
                    active_flag[params['parent']] = True

    return '\n'.join(ret)


def register():
    signals.get_generators.connect(initialize)
