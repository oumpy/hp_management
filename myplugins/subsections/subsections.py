"""
Subections Plugin for Pelican
-------

This plugin adds `subsections` attribute to article/page.
"""

from pelican import signals
from pelican.generators import ArticlesGenerator, StaticGenerator, PagesGenerator
import re
from bisect import bisect_right

import logging
logger = logging.getLogger(__name__)

class Section:
    def __init__(self, url, title, subsections = []):
        self.url, self.title = url, title
        if subsections == []:
            self.subsections = []
        else:
            self.subsections = subsections
section_id_fmt = 'pelican-subsections-{}'


def set_sections(obj):
    content = obj._content
    section_tags = re.finditer(r'(</|<)[hH][1-6](|\s+[^<>]*)>', content)
    comments = re.finditer(r'<!--.*?-->', content)
    # styles = re.finditer(r'<(style|STYLE)(|\s+[^<>]*)>.*?</(style|STYLE)>', content)
    comment_edges = []
    sectionlist = []
    idtag_stack = []
    for m in comments:
        comment_edges.append(m.start())
        comment_edges.append(m.end())
    n = 0
    contentlist = []
    prev_end = 0
    for m in section_tags:
        i = bisect_right(comment_edges, m.start())
        if i % 2: # in a comment
            continue
        else:
            tag = m.group()
            if tag[1] == '/':
                level = int(tag[3])
                idtag = idtag_stack.pop()
                title = content[prev_end:m.start()]
                sectionlist.append((level, idtag, title))
                contentlist.append(content[prev_end:m.start()])
                contentlist.append('</span>')
                contentlist.append(tag)
                prev_end = m.end()
            else:
                contentlist.append(content[prev_end:m.end()])
                idtag = section_id_fmt.format(n)
                n += 1
                idtag_stack.append(idtag)
                contentlist.append('<span id="{}">'.format(idtag))
                prev_end = m.end()
    else:
        contentlist.append(content[prev_end:])
    obj._content = ''.join(contentlist)

    setattr(obj, 'subsections', [])
    parents = [obj]
    for level, idtag, title in sectionlist:
        while len(parents) < level:
            parents.append(parents[-1])
        while len(parents) > level:
            parents.pop()
        section = Section(obj.url + '#' + idtag, title)
        parents[-1].subsections.append(section)
        parents.append(section)

def run_plugin(generators):
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            for article in generator.articles:
                set_sections(article)
        elif isinstance(generator, PagesGenerator):
            for page in generator.pages:
                set_sections(page)


def subsections_ancestorof(obj1, obj2):
    target_url = obj2.url
    pool = [obj1]
    while pool:
        x = pool.pop()
        if x.url == target_url:
            return True
        else:
            for s in x.subsections:
                pool.append(s)
    else:
        return False

def set_jinja(pelicanobj):
    pelicanobj.settings['JINJA_FILTERS']['subsections_ancestorof'] = subsections_ancestorof


def register():
    signals.initialized.connect(set_jinja)
    signals.all_generators_finalized.connect(run_plugin)
