"""
ArticleImageURL Plugin for Pelican
-------

This plugin applies Pelican URL convention to `image` & `social_image` metadata of articles.
"""

from __future__ import unicode_literals
from pelican import signals
import re

import logging

pelican_settings = {}
def _pelican_configure(pelicanobj):
    global pelican_settings
    for key in 'INTRASITE_LINK_REGEX', 'SITEURL':
        pelican_settings[key] = pelicanobj.settings[key]


def expand_link(link, article):
    link_regex = r"""^
        (?P<markup>)(?P<quote>)
        (?P<path>{0}(?P<value>.*))
        $""".format(article.settings['INTRASITE_LINK_REGEX'])
    links = re.compile(link_regex, re.X)
    return links.sub(
        lambda m: article._link_replacer(article.get_siteurl(), m),
        link)

def article_image_url(article):
    keys = ['image', 'social_image']
    for key in keys:
        if key in article.metadata.keys():
            article.metadata[key] = expand_link(article.metadata[key], article)
            logging.error(article.metadata[key])
            # logging.error(str(settings))
            # logging.error(settings['SITEURL'])
            # logging.error(str(article))
            logging.error(str(article._context['static_content']))
            #logging.error(str(generator._context['static_content']))
            #logging.error(str(generator.context['static_content']))
            # logging.error(str(article.metadata))
            # logging.error(str(article.static_content))
    if keys[0] in article.metadata.keys() and not keys[1] in article.metadata.keys():
        article.metadata[keys[1]] = article.metadata[keys[0]]
    elif (not keys[0] in article.metadata.keys()) and keys[1] in article.metadata.keys():
        article.metadata[keys[0]] = article.metadata[keys[1]]
    if keys[0] in article.metadata.keys():
        logging.error(article.metadata[keys[0]])

def run_plugin(generator, content):
    article_image_url(content)
    # for article in generator.articles:
    #     article_image_url(article)

def register_image_to_context():
    pass

def generate_context(static_generator, metadata):
    self.staticfiles = []
    linked_files = set(self.context['static_links'])
    found_files = self.get_files(self.settings['STATIC_PATHS'],
                                    exclude=self.settings['STATIC_EXCLUDES'],
                                    extensions=False)
    for f in linked_files | found_files:

        # skip content source files unless the user explicitly wants them
        if self.settings['STATIC_EXCLUDE_SOURCES']:
            if self._is_potential_source_path(f):
                continue

        static = self.readers.read_file(
            base_path=self.path, path=f, content_class=Static,
            fmt='static', context=self.context,
            preread_signal=signals.static_generator_preread,
            preread_sender=self,
            context_signal=signals.static_generator_context,
            context_sender=self)
        self.staticfiles.append(static)
        self.add_source_path(static, static=True)
    self._update_context(('staticfiles',))
    signals.static_generator_finalized.send(self)


def register():
    signals.initialized.connect(_pelican_configure)
    #signals.article_generator_finalized.connect(run_plugin)
    signals.static_generator_context(register_image_to_context)
    signals.article_generator_write_article.connect(run_plugin)
    #signals.content_object_init.connect(run_plugin)
