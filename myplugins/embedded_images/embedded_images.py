# -*- coding: utf-8 -*-
"""
EmbeddedImages Plugin for Pelican (in-house)
--------------------------------------------

Moves images embedded in article/page HTML as ``data:`` URIs -- the
output figures of Jupyter notebooks, Markdown-cell attachments, or any
base64 image pasted into a source -- out of the HTML and into files
under ``EMBEDDED_IMAGES_PATH`` in the output directory.  The HTML then
references them like any other static file, so pages shrink from
hundreds of KB to a few, browsers cache the figures, and feeds no
longer carry megabytes of base64.

How it works, and why no temporary files are involved:

- At ``content_object_init`` the plugin decodes every ``data:image/*``
  ``<img>`` in the raw content, keeps the bytes in memory, writes the
  file straight into the output directory (Pelican has already cleaned
  it by then), and rewrites the ``src`` to ``{static}/<path>/<hash>.<ext>``.
- A small stand-in object exposing ``url``/``save_as`` (what Pelican's
  link replacer needs) is registered in ``context['static_content']``
  under that path, so the ``{static}`` link is resolved by Pelican's
  standard machinery: relativized per page under ``RELATIVE_URLS``,
  absolute in feeds, correct under preview prefixes.  The object also
  carries ``data`` (the bytes) for other plugins, e.g. thumbnails.
- Pelican's static-link scan would otherwise try to read those paths
  from the content directory; the plugin removes them from
  ``context['static_links']`` before the StaticGenerator runs.

File names are content hashes, so identical images share one file and
URLs stay stable across edits.  When Pillow is available the ``<img>``
gets ``width``/``height`` attributes (prevents layout shift while
loading); ``loading="lazy"`` is added unless already present.

Settings:

- ``EMBEDDED_IMAGES_PATH`` : output subdirectory for the files
  (default ``'embedded_images'``).

Caveat: with ``LOAD_CONTENT_CACHE`` enabled, cached content skips
``content_object_init`` and its images would not be written; this site
does not use the cache.
"""

import base64
import hashlib
import os
import posixpath
import re
import io

from pelican import signals

import logging
logger = logging.getLogger(__name__)

try:
    from PIL import Image
except ImportError:  # pragma: no cover - optional
    Image = None

_IMG_TAG_RE = re.compile(r'<img\b[^>]*>', re.I | re.S)
_SRC_RE = re.compile(r'\bsrc\s*=\s*(["\'])(data:image/(?P<mime>[\w.+-]+)(?:;[\w.+-]+=[^;,"\']*)*;base64,(?P<data>[^"\']*))\1', re.I | re.S)
_ATTR_RE = re.compile(r'\b(width|height|loading)\s*=', re.I)

_MIME_EXT = {'jpeg': '.jpg', 'jpg': '.jpg', 'png': '.png', 'gif': '.gif',
             'webp': '.webp', 'svg+xml': '.svg', 'bmp': '.bmp'}


class GeneratedImage:
    """Stand-in for a static file that exists only in the output directory.

    Provides what Pelican's intrasite link replacer uses (``url``,
    ``save_as``, ``attach_to``) plus the image bytes for other plugins.
    """

    def __init__(self, url, data):
        self.url = url
        self.save_as = url
        self.data = data
        self.source_path = None

    def attach_to(self, content):
        # The location is fixed; {attach} behaves like {static}.
        pass


# path -> GeneratedImage, for the current build
_registry = {}


def _reset(pelican):
    _registry.clear()


def _dimensions(data):
    if Image is None:
        return None
    try:
        with Image.open(io.BytesIO(data)) as im:
            return im.size
    except Exception:
        return None


def _rewrite_tag(tag, instance, subdir, output_path):
    m = _SRC_RE.search(tag)
    if not m:
        return tag
    try:
        data = base64.b64decode(m.group('data'), validate=False)
    except (ValueError, TypeError):
        return tag
    if not data:
        return tag
    ext = _MIME_EXT.get(m.group('mime').lower(), '.bin')
    name = hashlib.sha1(data).hexdigest()[:16] + ext
    key = posixpath.join(subdir, name)

    obj = _registry.get(key)
    if obj is None:
        obj = GeneratedImage(key, data)
        _registry[key] = obj
        target = os.path.join(output_path, *key.split('/'))
        os.makedirs(os.path.dirname(target), exist_ok=True)
        if not os.path.exists(target) or os.path.getsize(target) != len(data):
            with open(target, 'wb') as f:
                f.write(data)
    instance._context['static_content'][key] = obj

    new_src = '{static}/' + key
    tag = tag[:m.start(2)] + new_src + tag[m.end(2):]

    extra = []
    if not _ATTR_RE.search(tag):
        size = _dimensions(data)
        if size:
            extra.append('width="{}" height="{}"'.format(*size))
        extra.append('loading="lazy"')
    if extra:
        end = tag.rstrip()
        if end.endswith('/>'):
            tag = end[:-2].rstrip() + ' ' + ' '.join(extra) + '/>'
        else:
            tag = end[:-1].rstrip() + ' ' + ' '.join(extra) + '>'
    return tag


def externalize(instance):
    content = getattr(instance, '_content', None)
    if not content or 'data:image/' not in content:
        return
    context = getattr(instance, '_context', None)
    if context is None or 'static_content' not in context:
        return
    settings = instance.settings
    subdir = settings.get('EMBEDDED_IMAGES_PATH', 'embedded_images').strip('/')
    output_path = settings['OUTPUT_PATH']
    instance._content = _IMG_TAG_RE.sub(
        lambda m: _rewrite_tag(m.group(0), instance, subdir, output_path),
        content)


def unlink_static_links(generator):
    """Keep the StaticGenerator from looking for our files in the content
    directory: it collects every {static} link of the articles and pages
    read so far.  Runs at page_generator_finalized, i.e. after articles
    and pages have been read and before the StaticGenerator runs."""
    links = generator.context.get('static_links')
    if links:
        links.difference_update(_registry.keys())


def register():
    signals.initialized.connect(_reset)
    signals.content_object_init.connect(externalize)
    signals.page_generator_finalized.connect(unlink_static_links)
