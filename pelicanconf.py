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
ARTICLE_SAVE_AS = ARTICLE_URL ='{category}/{date:%Y}/{date:%m}/{slug}.html'
PAGE_SAVE_AS = PAGE_URL ='{slug}.html'

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
PLUGINS = [
    'pelican-ipynb.markup',
    'render_math',
    'liquid_tags.youtube',
    'related_posts',
]

RELATED_POSTS_MAX = 3

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

CUSTOM_HTML_HEAD = "custom/html_head.html"
CUSTOM_HEADER_PAGE = "custom/header_page.html"
CUSTOM_HEADER_ARTICLE = "custom/header_article.html"
CUSTOM_ARTICLE_HEADERS = ("custom/article_header.html", "custom/open_in_colab_header.html",
                          "custom/toc_header.html", )
CUSTOM_INDEX_ARTICLE_HEADERS = ("custom/article_header.html",)

# Put taglist at end of articles, and use the default sharing button implementation.
CUSTOM_ARTICLE_FOOTERS = (
    "taglist.html", "sharing.html",
    "custom/utterances.html",
    "custom/related_posts.html",
)

CUSTOM_SCRIPTS_BASE = "custom/scripts_base.html"
CUSTOM_SCRIPTS_PAGE = "custom/page_showmodified_scripts.html"
CUSTOM_SCRIPTS_ARTICLE = "custom/scripts_article.html"
# SIDEBAR_HIDE_CATEGORIES = True

# Default sidebar template. Omit this setting for single column mode without sidebar.
SIDEBAR = "custom/sidebar.html"
CUSTOM_SIDEBAR_MIDDLES = ("custom/sb_links.html", "custom/sb_taglist.html", )
CUSTOM_SIDEBAR_BOTTOM = "custom/sb_twittertl.html"

TWITTER_TIMELINE_HEIGHT = 600
SIDEBAR_SIZE = 3

SOCIAL_SHARE_BUTTONS = (
    'twitter',
    'facebook',
    )

DISPLAY_RECENT_POSTS_ON_SIDEBAR=True

TWITTER_CARD = True
OPEN_GRAPH = True

# Read user's custom settings.
from content.contentconf import *

# Settings for Open Graph Properties
if not 'OPEN_GRAPH_ARTICLE_AUTHOR' in globals() and 'AUTHOR' in globals():
    OPEN_GRAPH_ARTICLE_AUTHOR = AUTHOR
if not 'DEFAULT_SOCIAL_IMAGE' in globals() and 'OPEN_GRAPH_IMAGE' in globals():
    DEFAULT_SOCIAL_IMAGE = OPEN_GRAPH_IMAGE
if not 'SOURCEREPOSITORY_URL' in globals():
    SOURCEREPOSITORY_URL = 'https://github.com/' + GITHUB_ACCOUNT + '/' + SOURCEREPOSITORY_NAME + '.git'
