# -*- coding: utf-8 -*- #
# give your GitHub Pages repository URL.
SITEREPOSITORY = 'https://github.com/oumpy/oumpy.github.io.git'

# If your site is available via HTTPS, make sure SITEURL begins with https://
SITEURL = 'oumpy.github.io'

# configuration for sitemap plugin
SITEMAP = {
    'format': 'xml',
    'priorities': {
        'articles': 0.5,
        'indexes': 0.5,
        'pages': 0.5
    },
    'changefreqs': {
        'articles': 'weekly',
        'indexes': 'weekly',
        'pages': 'weekly'
    }
}

# Following items are often useful when publishing
# Services
# GOOGLE_ANALYTICS = 'UA-12345678-9'
# DISQUS_SITENAME = 'johndoe'
