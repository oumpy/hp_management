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
    def __init__(self, url, title=None, subsections=[], active_pages=None, self_in_subsections=False):
        self.url = url
        self.title = title
        self.active_pages = active_pages
        self.self_in_subsections = self_in_subsections
        if subsections == []:
            self.subsections = []
        else:
            self.subsections = subsections
    AUTO = None
    DIVIDER = None

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
    dropdown_num = 0
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
            if node.title is None and node.url is not None and filter_url2obj:
                title = filter_url2obj(node.url).title
            else:
                title = node.title

            if s == FORWARD:
                if d >= depth:
                    subsections = []

                params['forward_line'] = len(ret)
                caret = ''
                params['/div'] = 0
                if node.url is None:
                    ret.append('<div class="dropdown-divider"></div>')
                    continue
                elif d == 0:
                    if subsections:
                        ret.append('<li class="nav-item dropdown dropdown-hover">')
                        if subsections and CARET:
                            caret = ' <span class="caret"></span>'
                        a_format = '<a class="nav-link dropdown-toggle" id="dropdown{}" href="{}" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">{}{}</a>'
                    else:
                        ret.append('<li class="nav-item">')
                        a_format = '<a class="nav-link" href="{1}">{2}{3}</a>'
                else:
                    if subsections:
                        ret.append('<div class="dropdown dropright">')
                        params['/div'] += 1
                        a_format = '<a class="dropdown-item dropdown-toggle" id="dropdown{}" href="{}" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">{}{}</a>'
                    else:
                        a_format = '<a class="dropdown-item" href="{1}">{2}{3}</a>'
                dropdown_num += 1
                ret.append(a_format.format(dropdown_num, '{}/{}'.format(rooturl, node.url), title, caret))

                if len(subsections) > 0:
                    ret.append('<div class="dropdown-menu" aria-labelledby="dropdown{}">'.format(dropdown_num))
                    params['/div'] += 1
                    pool.append((node, d, BACK, params))
                    for c in subsections[::-1]:
                        pool.append((c, d+1, FORWARD, {'parent':node.url}))
                    if node.self_in_subsections:
                        pool.append((MenuItem(MenuItem.DIVIDER), d+1, FORWARD, {'parent':node.url}))
                        pool.append((MenuItem(node.url,
                                              title=node.title,
                                              active_pages=node.active_pages,
                                              ),
                                     d+1, FORWARD, {'parent':node.url}))
                else:
                    pool.append((node, d, BACK, params))
            else: # s==BACK
                if params['/div'] > 0:
                    ret.append('</div>' * params['/div'])
                if d == 0:
                    ret.append('</li>')
                active_flag[node.url] |= bool(
                    hasattr(node,'active_pages') and node.active_pages and re.match(node.active_pages, page_url)
                )
                if active_flag[node.url]:
                    ret[params['forward_line']] = re.sub(r'(class="[^"]*)',r'\1 active', ret[params['forward_line']])
                    active_flag[params['parent']] = True

    return '\n'.join(ret)


def register():
    signals.get_generators.connect(initialize)
