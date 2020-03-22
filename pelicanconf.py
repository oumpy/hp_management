#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals
import datetime

AUTHOR = 'Python会'
SITENAME = '大阪大学医学部 Python会'
SITEURL = ''

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

# Blogroll
LINKS = (
        # ('Python会ブログ','https://oumedpython.hatenablog.com/'),
        #('大阪大学医学部', 'http://www.med.osaka-u.ac.jp/'),
        #('Python.org', 'http://python.org/'),
        #('Pelican（本サイトで使用）', 'http://getpelican.com/'),
        # ('Jinja2', 'http://jinja.pocoo.org/'),
        # ('You can modify those links in your config file', '#'),
        )

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
# THEME = './themes/voidy-bootstrap'
THEME = './theme/voidy-bootstrap'

# Theme customizations
# MINIMALXY_CUSTOM_CSS = 'theme/css/custom.css'
# MINIMALXY_FAVICON = 'favicon.ico'
# MINIMALXY_START_YEAR = 2018
# MINIMALXY_CURRENT_YEAR = datetime.date.today().year

# Author
AUTHOR_INTRO = u'大阪大学医学部所属のPython職人集団です'
AUTHOR_DESCRIPTION = u'Now is better than never'
# AUTHOR_AVATAR = '../images/'
# AUTHOR_WEB = 'https://twitter.com/oumed_python'

# Social
SOCIAL = (
    # ('facebook', ''),
    ('技術Blog (はてな)','https://oumedpython.hatenablog.com/'),
    ('Twitter', 'https://twitter.com/oumed_python'),
    ('E-mail', 'mailto:handai.python@gmail.com'),
    ('GitHub', 'https://github.com/oumpy'),
    ('YouTube', 'https://www.youtube.com/channel/UCh1eAeDCpsZeOh0Z9paNfHQ'),
)

###
### theme-specific settings
###

SITESUBTITLE ='Now is better than never.'
SITETAG = SITENAME

FONT_AWESOME_CDN_LINK = {
    'href': 'https://use.fontawesome.com/releases/v5.0.13/css/all.css',
    'integrity': 'sha384-DNOHZ68U8hZfKXOrtjWvjxusGo9WQnrNx2sqG0tfsghAvtVlRW3tvkXWZh58N9jp',
    'crossorigin': 'anonymous'
}

# Extra stylesheets, for bootstrap overrides or additional styling.
STYLESHEET_FILES = ("custom/pygment.css", "custom/voidybootstrap.css",)
CUSTOM_FOOTER = "custom/footer.html"
SKIP_COLOPHON = True

# Put taglist at end of articles, and use the default sharing button implementation.
CUSTOM_ARTICLE_FOOTERS = ("taglist.html", "sharing.html", )
CUSTOM_SCRIPTS_ARTICLE = "sharing_scripts.html"
# SIDEBAR_HIDE_CATEGORIES = True

# Settings for Twitter Timeline
CUSTOM_SIDEBAR_BOTTOM = "custom/sidebar_twittertimeline.html"
TWITTER_TIMELINE_URL = "https://twitter.com/oumed_python?ref_src=twsrc%5Etfw"
TWITTER_TIMELINE_HEIGHT = 720

# Default sidebar template. Omit this setting for single column mode without sidebar.
SIDEBAR = "custom/sidebar.html"
CUSTOM_SIDEBAR_MIDDLES = ("custom/sb_links.html", "custom/sb_taglist.html", )
SIDEBAR_SIZE = 3
SOCIAL_SHARE_BUTTONS = (
    'hatebu',
    'twitter',
    'facebook',
    'line',
    'pocket',
    'googleplus',
    )
TWITTER_USERNAME = 'oumed_python'

DISPLAY_RECENT_POSTS_ON_SIDEBAR=True

CUSTOM_SOCIAL_TITLE = "ソーシャル"
CUSTOM_CATEGORIES_TITLE = "記事カテゴリ"
CUSTOM_TAGS_TITLE = "タグ"
CUSTOM_LINKS_TITLE = "リンク"
