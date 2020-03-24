#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals
import datetime
import os
import sys
sys.path.append(os.curdir)

LOAD_CONTENT_CACHE = False

PATH = 'content'
STATIC_PATHS = ['images', 'extra', 'icon', 'data']
FAVICON = 'favicon.ico'
FAVICON_TYPE = 'image/vnd.microsoft.icon'
EXTRA_PATH_METADATA = {
    'icon/' + FAVICON: {'path': FAVICON},
}
OPEN_GRAPH_IMAGE = 'logo.jpg'

TIMEZONE = 'Asia/Tokyo'

DEFAULT_LANG = 'ja'
DATE_FORMATS = {
    'en': '%a, %d %b %Y',
    'ja': '%Y-%m-%d(%a)',
}

PAGE_ORDER_BY = 'page_order'

start_year = 2017
this_year = datetime.date.today().year

# ARTICLE_PATHS = ['articles']
ARTICLE_PATHS = [ 'articles/%dsy/' % y for y in range(start_year, this_year+2) ]
ARTICLE_SAVE_AS = 'articles/{date:%Y}/{date:%m}/{slug}.html'
ARTICLE_URL = 'articles/{date:%Y}/{date:%m}/{slug}.html'
PAGE_SAVE_AS = '{slug}.html'
PAGE_URL = '{slug}.html'

INDEX_SAVE_AS = 'articles.html'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

FEED_ALL_RSS = 'feeds/all.rss.xml'
FEED_ALL_ATOM = 'feeds/all.atom.xml'
RELATIVE_URLS = True

DEFAULT_PAGINATION = 5
PAGINATION_PATTERNS = (
    (1, '{url}', '{save_as}',),
    (2, '{base_name}/latests/{number}/', '{base_name}/latests/{number}/index.html'),
)
# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True
MARKUP = ('md', 'ipynb')

PLUGIN_PATHS = ['./plugins']
PLUGINS = ['pelican-ipynb.markup', 'render_math']

# if you create jupyter files in the content dir, snapshots are saved with the same
# metadata. These need to be ignored.
IGNORE_FILES = [".ipynb_checkpoints", '._*']

IPYNB_USE_METACELL = True
DISPLAY_PAGES_ON_MENU = True

USE_FOLDER_AS_CATEGORY = True

# Theme
THEME = './theme/voidy-bootstrap'

###
### theme-specific settings
###

FONT_AWESOME_CDN_LINK = {
    'href': 'https://use.fontawesome.com/releases/v5.0.13/css/all.css',
    'integrity': 'sha384-DNOHZ68U8hZfKXOrtjWvjxusGo9WQnrNx2sqG0tfsghAvtVlRW3tvkXWZh58N9jp',
    'crossorigin': 'anonymous'
}

# Extra stylesheets, for bootstrap overrides or additional styling.
STYLESHEET_FILES = ("pygment.css", "voidybootstrap.css", "voidybootstrap-custom.css")
CUSTOM_FOOTER = "custom/footer.html"
SKIP_COLOPHON = True

# Put taglist at end of articles, and use the default sharing button implementation.
CUSTOM_ARTICLE_FOOTERS = ("taglist.html", "sharing.html", )
CUSTOM_SCRIPTS_ARTICLE = "sharing_scripts.html"
# SIDEBAR_HIDE_CATEGORIES = True

# Default sidebar template. Omit this setting for single column mode without sidebar.
SIDEBAR = "custom/sidebar.html"
CUSTOM_SIDEBAR_MIDDLES = ("custom/sb_links.html", "custom/sb_taglist.html", )
CUSTOM_SIDEBAR_BOTTOM = "custom/sidebar_twittertimeline.html"

TWITTER_TIMELINE_HEIGHT = 600
SIDEBAR_SIZE = 3

SOCIAL_SHARE_BUTTONS = (
    'hatebu',
    'twitter',
    'facebook',
    'line',
    'pocket',
    'googleplus',
    )

DISPLAY_RECENT_POSTS_ON_SIDEBAR=True

# Read user's custom settings.
from content.contentconf import *
