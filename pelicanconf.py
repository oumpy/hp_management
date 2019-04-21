#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals
import datetime

AUTHOR = 'Python会'
SITENAME = '大阪大学医学部Python会'
SITEURL = ''

PATH = 'content'

TIMEZONE = 'Asia/Tokyo'

DEFAULT_LANG = 'ja'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Blogroll
LINKS = (('Pelican', 'http://getpelican.com/'),
         ('Python.org', 'http://python.org/'),
         ('Jinja2', 'http://jinja.pocoo.org/'),
         ('You can modify those links in your config file', '#'),)

# Social widget
SOCIAL = (('You can add links in your config file', '#'),
          ('Another social link', '#'),)

DEFAULT_PAGINATION = 10

# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True
MARKUP = ('md', 'ipynb')

PLUGIN_PATHS = ['./plugins']
PLUGINS = ['ipynb.markup']

# if you create jupyter files in the content dir, snapshots are saved with the same
# metadata. These need to be ignored.
IGNORE_FILES = [".ipynb_checkpoints", '._*']

IPYNB_USE_METACELL = True

# Theme
THEME = './themes/MinimalXY'

# Theme customizations
MINIMALXY_CUSTOM_CSS = 'theme/css/custom.css'
MINIMALXY_FAVICON = 'favicon.ico'
MINIMALXY_START_YEAR = 2018
MINIMALXY_CURRENT_YEAR = datetime.date.today().year

# Author
AUTHOR_INTRO = u'大阪大学医学部所属のPython職人集団です'
AUTHOR_DESCRIPTION = u'Now is better than never'
AUTHOR_AVATAR = '../images/'
AUTHOR_WEB = 'https://twitter.com/oumed_python'

# Services
GOOGLE_ANALYTICS = 'UA-12345678-9'
DISQUS_SITENAME = 'johndoe'

# Social
SOCIAL = (
    ('facebook', ''),
    ('twitter', 'https://twitter.com/oumed_python'),
    ('github', 'https://github.com/oumpy'),
    ('linkedin', 'https://www.youtube.com/channel/UCh1eAeDCpsZeOh0Z9paNfHQ'),
)
