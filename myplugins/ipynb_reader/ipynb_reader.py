# -*- coding: utf-8 -*-
"""
IPynb Reader Plugin for Pelican
-------------------------------

A lightweight, self-contained reader for Jupyter notebook (.ipynb) articles.

This plugin replaces the abandoned `pelican-jupyter` (markup mode), which
required nbconvert < 6 and pinned old versions of the whole Jupyter stack.
It depends only on packages Pelican itself already uses:

  - ``markdown`` : to render markdown cells (with the same extension set as
    regular .md articles, so MathJax handling via pelican-render-math,
    code highlighting etc. remain consistent site-wide),
  - ``pygments`` : to highlight code cells.

Compatibility with the previous system:

  - Metadata are read from a ``.nbdata`` file placed next to the notebook
    (same basename), written in the usual ``Key: Value`` format.
    The ``Subcells: [first, last]`` key is supported to select a range
    of cells, like pelican-jupyter's subcells feature.
  - The generated HTML mimics the CSS class structure of the old
    nbconvert 5 "basic" template (``cell``/``input_area``/``prompt``/
    ``output_subarea``/``highlight hl-ipython3`` ...), so the existing
    custom CSS and the "Open in Colab" button script keep working without
    modification.
"""

import ast
import base64
import json
import os
import re

from pelican import signals
from pelican.readers import BaseReader
from pelican.utils import pelican_open

from markdown import Markdown
from pygments import highlight as pygments_highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import TextLexer, get_lexer_by_name
from pygments.util import ClassNotFound

import logging
logger = logging.getLogger(__name__)


ANSI_RE = re.compile(r'\x1b\[[0-9;]*[a-zA-Z]')

# Priority order of MIME types for rich outputs.
DISPLAY_DATA_PRIORITY = [
    'text/html',
    'image/svg+xml',
    'image/png',
    'image/jpeg',
    'text/latex',
    'text/markdown',
    'text/plain',
]


def strip_ansi(text):
    return ANSI_RE.sub('', text)


def escape_html(text):
    return (
        text.replace('&', '&amp;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
    )


def join_source(source):
    if isinstance(source, list):
        return ''.join(source)
    return source or ''


class IPynbHtmlBuilder(object):
    """Convert a notebook (nbformat 4 dict) to an HTML fragment."""

    def __init__(self, settings, markdown_settings):
        self.settings = settings
        self._md = Markdown(**markdown_settings)

    # -- cell renderers ----------------------------------------------------

    def render_markdown(self, text, attachments=None):
        if attachments:
            for name, data in attachments.items():
                for mime, b64 in data.items():
                    text = text.replace(
                        'attachment:' + name,
                        'data:{};base64,{}'.format(mime, b64))
        self._md.reset()
        return self._md.convert(text)

    def highlight_code(self, source, language='ipython3'):
        try:
            if language in ('ipython3', 'ipython', 'python', 'python3', None, ''):
                lexer = get_lexer_by_name('python3')
            else:
                lexer = get_lexer_by_name(language)
        except ClassNotFound:
            lexer = TextLexer()
        formatter = HtmlFormatter(cssclass='highlight hl-ipython3')
        return pygments_highlight(source, lexer, formatter)

    def markdown_cell(self, cell):
        html = self.render_markdown(join_source(cell.get('source')),
                                    cell.get('attachments'))
        return (
            '<div class="cell border-box-sizing text_cell rendered">'
            '<div class="prompt input_prompt"></div>'
            '<div class="inner_cell">'
            '<div class="text_cell_render border-box-sizing rendered_html">'
            '{}'
            '</div></div></div>\n'.format(html)
        )

    def raw_cell(self, cell):
        # Raw cells are passed through untouched (same as nbconvert with
        # raw_mimetype unset rendering to nothing useful; we escape them).
        return '<pre class="raw_cell">{}</pre>\n'.format(
            escape_html(join_source(cell.get('source'))))

    def code_cell(self, cell, language):
        count = cell.get('execution_count')
        prompt_in = 'In [{}]:'.format(count if count is not None else ' ')
        parts = [
            '<div class="cell border-box-sizing code_cell rendered">',
            '<div class="input">',
            '<div class="prompt input_prompt">{}</div>'.format(prompt_in),
            '<div class="inner_cell"><div class="input_area">',
            self.highlight_code(join_source(cell.get('source')), language),
            '</div></div>',
            '</div>',  # /input
        ]
        outputs = cell.get('outputs') or []
        if outputs:
            parts.append('<div class="output_wrapper"><div class="output">')
            for output in outputs:
                parts.append(self.render_output(output))
            parts.append('</div></div>')
        parts.append('</div>\n')
        return ''.join(parts)

    # -- output renderers --------------------------------------------------

    def render_output(self, output):
        otype = output.get('output_type')
        if otype == 'stream':
            css = ('output_stdout' if output.get('name') != 'stderr'
                   else 'output_stderr')
            return (
                '<div class="output_area">'
                '<div class="prompt"></div>'
                '<div class="output_subarea output_text output_stream {}">'
                '<pre>{}</pre></div></div>'.format(
                    css, escape_html(strip_ansi(join_source(output.get('text')))))
            )
        elif otype == 'error':
            tb = '\n'.join(output.get('traceback') or [])
            return (
                '<div class="output_area">'
                '<div class="prompt"></div>'
                '<div class="output_subarea output_text output_error">'
                '<pre>{}</pre></div></div>'.format(
                    escape_html(strip_ansi(tb)))
            )
        elif otype in ('execute_result', 'display_data'):
            count = output.get('execution_count')
            if otype == 'execute_result' and count is not None:
                prompt = 'Out[{}]:'.format(count)
            else:
                prompt = ''
            body = self.render_data(output.get('data') or {})
            return (
                '<div class="output_area">'
                '<div class="prompt output_prompt">{}</div>'
                '{}</div>'.format(prompt, body)
            )
        return ''

    def render_data(self, data):
        for mime in DISPLAY_DATA_PRIORITY:
            if mime not in data:
                continue
            content = data[mime]
            if isinstance(content, list):
                content = ''.join(content)
            if mime == 'text/html':
                return ('<div class="output_subarea output_html rendered_html">'
                        '{}</div>'.format(content))
            elif mime == 'image/svg+xml':
                return ('<div class="output_subarea output_svg">'
                        '{}</div>'.format(content))
            elif mime in ('image/png', 'image/jpeg'):
                b64 = content.replace('\n', '')
                return ('<div class="output_subarea output_png">'
                        '<img src="data:{};base64,{}" alt=""/></div>'.format(mime, b64))
            elif mime == 'text/latex':
                # Left as-is for MathJax to typeset.
                return ('<div class="output_subarea output_latex math">'
                        '{}</div>'.format(content))
            elif mime == 'text/markdown':
                return ('<div class="output_subarea output_markdown rendered_html">'
                        '{}</div>'.format(self.render_markdown(content)))
            elif mime == 'text/plain':
                return ('<div class="output_subarea output_text">'
                        '<pre>{}</pre></div>'.format(escape_html(content)))
        return '<div class="output_subarea"></div>'

    # -- whole notebook ----------------------------------------------------

    def build(self, notebook, subcells=None):
        if notebook.get('nbformat', 0) < 4:
            raise ValueError('Only nbformat >= 4 notebooks are supported.')
        metadata = notebook.get('metadata') or {}
        language = (metadata.get('kernelspec') or {}).get('language') \
            or (metadata.get('language_info') or {}).get('name') or 'python'
        cells = notebook.get('cells') or []
        if subcells is not None:
            cells = cells[slice(*subcells)]
        parts = []
        for cell in cells:
            ctype = cell.get('cell_type')
            if ctype == 'markdown':
                parts.append(self.markdown_cell(cell))
            elif ctype == 'code':
                parts.append(self.code_cell(cell, language))
            elif ctype == 'raw':
                parts.append(self.raw_cell(cell))
        return '<div class="jupyter-notebook">\n{}</div>'.format(''.join(parts))


class IPynbReader(BaseReader):
    """Reader for .ipynb files with .nbdata metadata files."""

    enabled = True
    file_extensions = ['ipynb']

    def read(self, source_path):
        # --- metadata from the companion .nbdata file
        nbdata_path = os.path.splitext(source_path)[0] + '.nbdata'
        if not os.path.exists(nbdata_path):
            raise Exception(
                'Missing metadata file {} for notebook {}.'.format(
                    nbdata_path, source_path))
        metadata = {}
        subcells = None
        with pelican_open(nbdata_path) as text:
            for line in text.splitlines():
                if not line.strip():
                    continue
                if ':' not in line:
                    logger.warning('Unrecognized nbdata line in %s: %r',
                                   nbdata_path, line)
                    continue
                key, value = line.split(':', 1)
                key = key.strip().lower()
                value = value.strip()
                if key == 'subcells':
                    subcells = ast.literal_eval(value)
                else:
                    metadata[key] = self.process_metadata(key, value)

        # --- notebook body
        with pelican_open(source_path) as text:
            notebook = json.loads(text)

        markdown_settings = dict(self.settings.get('MARKDOWN') or {})
        markdown_settings.setdefault('output_format', 'html5')
        builder = IPynbHtmlBuilder(self.settings, markdown_settings)
        content = builder.build(notebook, subcells)
        return content, metadata


def add_reader(readers):
    readers.reader_classes['ipynb'] = IPynbReader


def register():
    signals.readers_init.connect(add_reader)
