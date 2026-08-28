"""
MakeMenu Plugin for Pelican
-------

This plugin provides the `resolve_menu` Jinja2 filter, which resolves the
site's `ADD_ON_MENU` declaration (a tree of `MenuItem` objects, page/article
objects, or URL strings) into a plain tree of dicts:

    {'url': ..., 'title': ..., 'active': bool, 'divider': bool,
     'children': [...]}

The HTML markup itself is rendered by the theme's Jinja2 macros
(`templates/includes/menu.html`); this module contains no HTML.
Relative URL handling is left to Pelican: templates simply prefix
`{{ SITEURL }}`, which Pelican relativizes per page when
`RELATIVE_URLS = True`.

The traversal is an iterative depth-first search with an explicit stack
and FORWARD/BACK phases, preserved from the original implementation:
each node is visited twice, FORWARD to construct its entry (and push its
subsections), BACK to settle its `active` state after all descendants
have been processed, propagating activity to the parent via
`active_flag`.

Notes:
- `subsections=MenuItem.AUTO` (None) expands the subsections attached by
  the `subsections` plugin, looked up via the `url2obj` filter provided by
  the `path2obj` plugin.
- `active` is true for the page itself, for nodes whose `active_pages`
  regex matches the current page URL, and propagates from children to
  their ancestors.
"""

from __future__ import unicode_literals
from pelican import signals
from collections import defaultdict
import re

import logging
logger = logging.getLogger(__name__)

filter_url2obj = None

def initialize(pelicanobj):
    global filter_url2obj
    settings = pelicanobj.settings
    if 'url2obj' in settings['JINJA_FILTERS'].keys():
        filter_url2obj = settings['JINJA_FILTERS']['url2obj']
    settings['JINJA_FILTERS']['resolve_menu'] = resolve_menu


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


def _divider():
    return {'url': None, 'title': None, 'active': False,
            'divider': True, 'children': []}


def resolve_menu(add_on_menu, page_url, depth=1):
    """Resolve ADD_ON_MENU into a plain tree of dicts for the template.

    add_on_menu : list/tuple of MenuItem objects (or URL strings)
    page_url    : URL of the page being rendered (for active detection)
    depth       : menu hierarchy depth (MENU_STEPS)
    """
    roots = []
    FORWARD, BACK = 0, 1
    for obj in add_on_menu:
        pool = [(obj, 0, FORWARD, {'parent': None, 'siblings': roots})]
        active_flag = defaultdict(bool)
        active_flag[page_url] = True
        while pool:
            node, d, s, params = pool.pop()
            if isinstance(node, str):
                if filter_url2obj:
                    node = filter_url2obj(node)
                else:
                    logger.error('You need path2obj for \'{}\' in submenu.'.format(node))
                    continue
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
                if node.url is None:
                    params['siblings'].append(_divider())
                    continue
                if d >= depth:
                    subsections = []

                item = {'url': node.url, 'title': title, 'active': False,
                        'divider': False, 'children': []}
                params['siblings'].append(item)
                params['item'] = item
                pool.append((node, d, BACK, params))

                if len(subsections) > 0:
                    for c in subsections[::-1]:
                        pool.append((c, d + 1, FORWARD,
                                     {'parent': node.url,
                                      'siblings': item['children']}))
                    if node.self_in_subsections:
                        pool.append((MenuItem(MenuItem.DIVIDER), d + 1, FORWARD,
                                     {'parent': node.url,
                                      'siblings': item['children']}))
                        pool.append((MenuItem(node.url,
                                              title=node.title,
                                              active_pages=node.active_pages,
                                              ),
                                     d + 1, FORWARD,
                                     {'parent': node.url,
                                      'siblings': item['children']}))
            else:  # s == BACK
                active_flag[node.url] |= bool(
                    hasattr(node, 'active_pages') and node.active_pages and re.match(node.active_pages, page_url)
                )
                if active_flag[node.url]:
                    params['item']['active'] = True
                    active_flag[params['parent']] = True

    return roots


def register():
    signals.get_generators.connect(initialize)
