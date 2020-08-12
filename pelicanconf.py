#!/usr/bin/env python
# -*- coding: utf-8 -*- #
import datetime
import copy
import os
import sys
sys.path.append(os.curdir)

BOOTSTRAP_VERSION = '4.5.0'

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

# PAGE_ORDER_BY = 'page_order'

ARTICLE_PATHS = ['articles']
ARTICLE_SAVE_AS = ARTICLE_URL ='{category}/{date:%Y}/{date:%m}/{slug}.html'
PAGE_SAVE_AS = PAGE_URL ='{slug}.html'
CATEGORY_SAVE_AS = CATEGORY_URL = '{slug}.html'
ARTICLE_EXCLUDES_DIRNAMES = PAGE_EXCLUDES_DIRNAMES = ['attach', 'images']
INDEX_SAVE_AS = 'articles.html'

SLUGIFY_SOURCE = 'basename'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

FEED_ALL_RSS = 'feeds/all.rss.xml'
FEED_ALL_ATOM = 'feeds/all.atom.xml'
RELATIVE_URLS = True

TAGS_URL = TAGS_SAVE_AS = 'tags.html'
AUTHORS_URL = AUTHORS_SAVE_AS = 'authors.html'

DEFAULT_PAGINATION = 10
PAGINATION_PATTERNS = (
    (1, '{url}', '{save_as}',),
    (2, '{base_name}/latests/{number}/', '{base_name}/latests/{number}/index.html'),
)
CUSTOM_CONTENT_TOP_CATEGORY = "custom/content_top_category.html"

# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True
MARKUP = ['md', 'ipynb']

PLUGIN_PATHS = ['./plugins', './myplugins']
PLUGINS = [
    'pelican-ipynb.markup',
    'render_math',
    'tag_cloud',
    'related_posts',
    'autosummary', 'summary', # this order is important!
    'category_names',
    'shortcodes',
    'apply_jinja2',
    'path2obj',
    'subsections',
    'makemenu',
    'pelican-sass',
    'excludes_dirnames',
    'skiptags',
]

RELATED_POSTS_MAX = 3

SHORTCODES = {
    'youtube': '''\
<p><span class="videobox">
  <iframe width="{{width|default(640)}}" height="{{height|default(390)}}"
    src="https://www.youtube.com/embed/{{id}}"
    frameborder="0" webkitAllowFullScreen mozallowfullscreen allowFullScreen>
  </iframe></span></p>''',
    'embed': '''\
<p><span class="videobox">
  <iframe width={{width|default(640)}}" height="{{height|default(390)}}"
    src="{{src}}"
    frameborder="0" webkitAllowFullScreen mozallowfullscreen allowFullScreen>
  </iframe></span></p>''',
}

# if you create jupyter files in the content dir, snapshots are saved with the same
# metadata. These need to be ignored.
IGNORE_FILES = [".ipynb_checkpoints", '._*']

# IPYNB_USE_METACELL = True
# DISPLAY_PAGES_ON_MENU = True
# USE_FOLDER_AS_CATEGORY = True

# Theme
THEME = './theme/voidy-bootstrap'

###
### theme-specific settings
###

FONT_AWESOME_LINK = {
    'href': 'https://use.fontawesome.com/releases/v5.13.0/css/all.css',
    'integrity': 'sha384-Bfad6CLCknfcloXFOyFnlgtENryhrpZCe29RTifKEixXQZ38WheV+i/6YWSzkz3V',
    'crossorigin': 'anonymous'
}

# Extra stylesheets, for bootstrap overrides or additional styling.
STYLESHEET_FILES = [
    "pygment.css",
    # "voidybootstrap.css",
    "theme.css",
    "voidybootstrap-custom.css",
]
CUSTOM_FOOTER = "custom/footer.html"
SKIP_COLOPHON = True

CUSTOM_HTML_HEAD = "custom/html_head.html"
CUSTOM_HEADER_PAGE = "custom/header_page.html"
CUSTOM_HEADER_ARTICLE = "custom/header_article.html"
CUSTOM_ARTICLE_HEADERS = [
    "custom/article_header.html",
    "custom/open_in_colab_header.html",
    "custom/toc_header.html",
]
CUSTOM_INDEX_ARTICLE_HEADERS = [
    "custom/article_header.html",
]

# Put taglist at end of articles, and use the default sharing button implementation.
CUSTOM_ARTICLE_FOOTERS = [
    "taglist.html",
    "sharing.html",
    "custom/utterances.html",
    "custom/related_posts.html",
]

CUSTOM_SCRIPTS_BASE = "custom/scripts_base.html"
CUSTOM_SCRIPTS_PAGE = "custom/page_showmodified_scripts.html"
CUSTOM_SCRIPTS_ARTICLE = "custom/scripts_article.html"
# SIDEBAR_HIDE_CATEGORIES = True

# Default sidebar template. Omit this setting for single column mode without sidebar.
SIDEBAR = "custom/sidebar.html"
CUSTOM_SIDEBAR_MIDDLES = [
    "custom/sb_google_cse.html",
    "custom/sb_social.html",
    "custom/sb_recentposts.html",
    "custom/sb_tagcloud.html",
    "custom/sb_support.html"
]
CUSTOM_SIDEBAR_BOTTOM = "custom/sb_twittertl.html"
SIDEBAR_HIDE_FEEDS = True

CUSTOM_CONTENT_TOP_ARCHIVES = "custom/content_top_archives.html"

TWITTER_TIMELINE_HEIGHT = 600
SIDEBAR_SIZE = 3

SOCIAL_SHARE_BUTTONS = [
    'twitter',
    'facebook',
]

DISPLAY_RECENT_POSTS_ON_SIDEBAR=True

TWITTER_CARD = True
OPEN_GRAPH = True


# Read user's custom settings.
import tools.lib.pelicanns
globals_copy = copy.copy(globals())
for k, v in globals_copy.items():
    setattr(tools.lib.pelicanns, k, v)
from content.contentconf import *

# Settings for Open Graph Properties
if not 'OPEN_GRAPH_ARTICLE_AUTHOR' in globals() and 'AUTHOR' in globals():
    OPEN_GRAPH_ARTICLE_AUTHOR = AUTHOR
if not 'DEFAULT_SOCIAL_IMAGE' in globals() and 'OPEN_GRAPH_IMAGE' in globals():
    DEFAULT_SOCIAL_IMAGE = OPEN_GRAPH_IMAGE
if not 'SOURCEREPOSITORY_URL' in globals():
    SOURCEREPOSITORY_URL = 'https://github.com/' + GITHUB_ACCOUNT + '/' + SOURCEREPOSITORY_NAME + '.git'
