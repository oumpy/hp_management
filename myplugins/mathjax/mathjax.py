# -*- coding: utf-8 -*-
"""
MathJax Plugin for Pelican (in-house, MathJax v4)
-------------------------------------------------

Renders LaTeX math in Markdown content (including the markdown cells of
.ipynb articles processed by myplugins/ipynb_reader) with MathJax v4.

This replaces the pelican-render-math plugin, which hardcodes a MathJax
2.7.3 loader (2017) whose configuration mechanism (MathJax.Hub.Config)
is incompatible with MathJax v3/v4.  The valuable part of render-math is
its Markdown extension protecting math from the Markdown processor;
that extension is vendored here unchanged
(pelican_mathjax_markdown_extension.py, AGPL -- see LICENSE), while the
loader script is rewritten for the current MathJax component system.

How it works:

- The Markdown extension wraps ``$...$`` (inline) as
  ``<span class="math">\\(...\\)</span>`` and ``$$...$$`` /
  ``\\begin{env}...`` (display) as ``<div class="math">...</div>``,
  and appends the loader script to any document containing math.
- The loader sets ``window.MathJax`` (TeX input configuration) and then
  loads the ``tex-chtml`` component from jsDelivr.  It is idempotent, so
  documents assembled from many fragments (e.g. notebook cells) load
  MathJax only once.
- Summaries that contain math get the loader appended too (and a
  truncated last formula is completed), so listing pages render math
  correctly -- same behavior as render-math.

Settings:

- ``MATHJAX_SOURCE`` (optional): URL of the MathJax component to load.
  Default: ``https://cdn.jsdelivr.net/npm/mathjax@4/tex-chtml.js``
  (major-pinned; patch releases are picked up automatically).
"""

import functools

from bs4 import BeautifulSoup

from pelican import generators, signals

from .pelican_mathjax_markdown_extension import PelicanMathJaxExtension

import logging
logger = logging.getLogger(__name__)

DEFAULT_MATHJAX_SOURCE = 'https://cdn.jsdelivr.net/npm/mathjax@4/tex-chtml.js'

# Raw JavaScript (the markdown extension wraps it in a <script> element).
# Idempotent: only the first fragment on a page actually loads MathJax.
MATHJAX_SCRIPT_TEMPLATE = """
if (!window._pelican_mathjax_loading) {{
    window._pelican_mathjax_loading = true;
    window.MathJax = Object.assign({{
        tex: {{
            inlineMath: [['\\\\(', '\\\\)']],
            displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
            processEscapes: true
        }}
    }}, window.MathJax || {{}});
    var mathjaxscript = document.createElement('script');
    mathjaxscript.id = 'mathjaxscript_pelican';
    mathjaxscript.defer = true;
    mathjaxscript.src = '{source}';
    document.head.appendChild(mathjaxscript);
}}
"""

mathjax_script = None


def pelican_init(pelicanobj):
    global mathjax_script
    settings = pelicanobj.settings
    source = settings.get('MATHJAX_SOURCE', DEFAULT_MATHJAX_SOURCE)
    mathjax_script = MATHJAX_SCRIPT_TEMPLATE.format(source=source)

    config = {
        'mathjax_script': mathjax_script,
        'math_tag_class': 'math',
        'auto_insert': True,
    }
    settings['MARKDOWN'].setdefault('extensions', []).append(
        PelicanMathJaxExtension(config))


_MATH_DELIMITERS = ('\\(', '\\[', '$$')


def process_summary(article):
    """Complete a truncated last formula and add the MathJax loader
    to summaries containing math (for index/listing pages).

    Must run after any plugin that rewrites the summary (autosummary):
    that one strips the ``<span class="math">`` wrappers but leaves the
    TeX delimiters in the text, which MathJax typesets on its own as
    long as the loader is present."""
    summary = article.summary
    if not summary:
        return
    summary_parsed = BeautifulSoup(summary, 'html.parser')
    math = summary_parsed.find_all(class_='math')
    if not math and not any(d in summary_parsed.get_text()
                            for d in _MATH_DELIMITERS):
        return

    if math:
        last_math_text = math[-1].get_text()
        if len(last_math_text) > 3 and last_math_text[-3:] == '...':
            content_parsed = BeautifulSoup(article._content, 'html.parser')
            full_text = content_parsed.find_all(class_='math')[len(math) - 1].get_text()
            math[-1].string = '%s ...' % full_text
            summary = summary_parsed.decode()

    # clear memoization cache
    if isinstance(article.get_summary, functools.partial):
        memoize_instance = article.get_summary.func.__self__
        memoize_instance.cache.clear()

    article.metadata['summary'] = (
        f"{summary}<script type='text/javascript'>{mathjax_script}</script>"
    )


def process_summaries(content_generators):
    for generator in content_generators:
        if isinstance(generator, generators.ArticlesGenerator):
            for article in (
                generator.articles + generator.translations + generator.drafts
            ):
                process_summary(article)


def register():
    signals.initialized.connect(pelican_init)
    signals.all_generators_finalized.connect(process_summaries)
