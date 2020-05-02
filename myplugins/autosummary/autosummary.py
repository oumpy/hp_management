"""
AutoSummary
-------

This plugin generates summary from article, with customizable tag-selections.
"""

from __future__ import unicode_literals
from bs4 import BeautifulSoup
from pelican import signals
from pelican.generators import ArticlesGenerator, StaticGenerator, PagesGenerator
import re

def initialized(pelican):
    from pelican.settings import DEFAULT_CONFIG
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_KEEP_HTMLTAGS',
                              ('a', 'font', 's', 'strong', 'em', 'u', 'b'))
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_REPLACE_HTMLTAGS',
                              {'h[1-9]' : 'strong'})
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_MAX_LENGTH', 140)
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_MIN_LENGTH', 1)
    if pelican:
        pelican.settings.setdefault('AUTOSUMMARY_KEEP_HTMLTAGS',
                                    ('a', 'font', 's', 'strong', 'em', 'u', 'b'))
        pelican.settings.setdefault('AUTOSUMMARY_REPLACE_HTMLTAGS',
                                    {'h[1-9]' : 'strong'})
        pelican.settings.setdefault('AUTOSUMMARY_MAX_LENGTH', 140)
        pelican.settings.setdefault('AUTOSUMMARY_MIN_LENGTH', 1)


def make_autosummary(content, settings):
    content = str(BeautifulSoup(content, 'html.parser'))
    content = re.compile(r'([\s　¶]|&#182;)+').sub(' ', content) # &#182; appears in converted jupyter notebooks.
    content = re.compile(r'<(style|STYLE)(|\s+[^<>]*)>.*?</(style|STYLE)>').sub(' ', content) # for jupyter. in general, probably buggy.
    content = re.compile(r'[\s　]+').sub(' ', content)
    cur_pos = 0
    length = 0
    max_length = settings['AUTOSUMMARY_MAX_LENGTH']
    htmltags = settings['AUTOSUMMARY_KEEP_HTMLTAGS']
    replacetags = settings['AUTOSUMMARY_REPLACE_HTMLTAGS']
    tag_stack = []
    summary = ''
    for m in re.finditer(r'<[^>]*?>', content + '<>'): # HTML tags
        start = m.start()
        end = m.end()
        tag = m.group()
        section_len = start - cur_pos - content.count(' ', cur_pos, start)
        if length + section_len <= max_length:
            tag_sp = re.match(r'(</|<)([^\s/>]*)', tag).group(2).lower()
            for tag_regex in htmltags:
                if re.fullmatch(tag_regex, tag_sp):
                    if tag[1] == '/':
                        tag_stack.pop()
                    else:
                        tag_stack.append(tag_sp)
                    summary += content[cur_pos: end]
                    break
            else:
                for reptag_regex, newtag in replacetags.items():
                    if re.fullmatch(reptag_regex, tag_sp):
                        summary += content[cur_pos: start]
                        if tag[1] == '/':
                            tag_stack.pop()
                            summary += '</%s>' % newtag
                        else:
                            tag_stack.append(newtag)
                            summary += '<%s>' % newtag
                        break
                else:
                    summary += content[cur_pos: start]
            length += section_len
            cur_pos = end
        else:
            w = max_length - length
            i = 0
            p = cur_pos
            while (i < w):
                if content[p] != ' ':
                    i += 1
                p += 1
            summary += content[cur_pos: p]
            for tag_sp in tag_stack[::-1]:
                summary += '</%s>' % tag_sp
            summary += '....'
            break
    return summary


def contentlen(content):
    content = str(BeautifulSoup(content, 'html.parser'))
    content = re.compile(r'([\s　¶]|&#182;)+').sub('', content) # &#182; appears in converted jupyter notebooks.
    content = re.compile(r'<(style|STYLE)(|\s+[^<>]*)>.*?</(style|STYLE)>').sub('', content) # for jupyter. in general, probably buggy.
    return len(re.compile(r'<[^>]*?>').sub('', content))


def extract_summary(instance):
    min_length = instance.settings['AUTOSUMMARY_MIN_LENGTH']
    # if summary is already specified, use it
    # if there is no content, there's nothing to do
    if ('summary' in instance.metadata and
        (contentlen(instance.metadata['summary']) >= min_length or ((not instance._content) and not hasattr(instance, '_summary')))):
        instance.has_summary = True
        return
    elif hasattr(instance, '_summary') and (contentlen(instance._summary) >= min_length or (not instance._content)):
        pre_summary = instance._summary
    elif not instance._content:
        instance.has_summary = False
        return
    else:
        pre_summary = instance._content

    pre_summary = instance._update_content(pre_summary, instance.settings['SITEURL'])
    summary = make_autosummary(pre_summary, instance.settings)

    # default_status was added to Pelican Content objects after 3.7.1.
    # Its use here is strictly to decide on how to set the summary.
    # There's probably a better way to do this but I couldn't find it.
    if hasattr(instance, 'default_status'):
        instance.metadata['summary'] = summary
    else:
        instance._summary = summary
    instance.has_summary = True


def run_plugin(generators):
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            for article in generator.articles:
                extract_summary(article)
        elif isinstance(generator, PagesGenerator):
            for page in generator.pages:
                extract_summary(page)


def register():
    signals.initialized.connect(initialized)
    try:
        signals.all_generators_finalized.connect(run_plugin)
    except AttributeError:
        # NOTE: This results in #314 so shouldn't really be relied on
        # https://github.com/getpelican/pelican-plugins/issues/314
        signals.content_object_init.connect(extract_summary)
