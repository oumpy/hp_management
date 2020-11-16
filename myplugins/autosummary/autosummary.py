"""
AutoSummary
-------

This plugin generates summary from article, with customizable tag-selections.
"""

from bs4 import BeautifulSoup
from pelican import signals
from pelican.generators import ArticlesGenerator, StaticGenerator, PagesGenerator
import re

import logging
logger = logging.getLogger(__name__)

def initialized(pelican):
    from pelican.settings import DEFAULT_CONFIG
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_KEEP_HTMLTAGS',
                              ('a', 'font', 's', 'strong', 'em', 'u', 'b'))
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_REPLACE_HTMLTAGS',
                              {'h[2-9]' : 'strong'})
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_DELETE_HTMLTAGS', ['h1'])
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_MAX_LENGTH', 140)
    DEFAULT_CONFIG.setdefault('AUTOSUMMARY_MIN_LENGTH', 1)
    if pelican:
        pelican.settings.setdefault('AUTOSUMMARY_KEEP_HTMLTAGS',
                                    ('a', 'font', 's', 'strong', 'em', 'u', 'b'))
        pelican.settings.setdefault('AUTOSUMMARY_REPLACE_HTMLTAGS',
                                    {'h[2-9]' : 'strong'})
        pelican.settings.setdefault('AUTOSUMMARY_DELETE_HTMLTAGS', ['h1'])
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
    deletetags = settings['AUTOSUMMARY_DELETE_HTMLTAGS']
    tag_stack = []
    summary = ''
    comment_level = 0
    for m in re.finditer(r'<[^>]*?>', content + '<>'): # HTML tags
        start = m.start()
        end = m.end()
        tag = m.group()
        section_len = start - cur_pos - content.count(' ', cur_pos, start)
        if comment_level > 0 or length + section_len <= max_length:
            tag_sp = re.match(r'(</|<)([^\s/>]*)', tag).group(2).lower()
            if not comment_level:
                length += section_len
            for tag_regex in htmltags:
                if re.fullmatch(tag_regex, tag_sp):
                    if tag[1] == '/':
                        tag_stack.pop()
                    else:
                        tag_stack.append(tag_sp)
                    if not comment_level:
                        summary += content[cur_pos: end]
                    break
            else:
                for reptag_regex, newtag in replacetags.items():
                    if re.fullmatch(reptag_regex, tag_sp):
                        if not comment_level:
                            summary += content[cur_pos: start]
                        if tag[1] == '/':
                            tag_stack.pop()
                            if not comment_level:
                                summary += '</%s>' % newtag
                        else:
                            tag_stack.append(newtag)
                            if not comment_level:
                                summary += '<%s>' % newtag
                        break
                else:
                    for deltag_regex in deletetags:
                        if re.fullmatch(deltag_regex, tag_sp):
                            if tag[1] == '/':
                                tag_stack.pop()
                                comment_level -= 1
                            else:
                                tag_stack.append(tag_sp)
                                comment_level += 1
                            break
                    else:
                        if not comment_level:
                            summary += content[cur_pos: start]
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
                for deltag_regex in deletetags:
                    if re.fullmatch(deltag_regex, tag_sp):
                        comment_level -= 1
                        break
                else:
                    if not comment_level:
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
    if hasattr(instance, '_summary') and (contentlen(instance._summary) >= min_length or (not instance._content)):
        pre_summary = instance._summary
    elif not instance._content:
    # if there is no content, there's nothing to do
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
