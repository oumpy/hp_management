# -*- coding: utf-8 -*-
"""
Thumbnails Plugin for Pelican (in-house)
----------------------------------------

Gives every article a thumbnail image for use in article lists (index,
category, tag and author pages), in the spirit of WordPress "featured
images".  The image is chosen, in order of preference, from:

1. the article's ``thumbnail`` metadata, written exactly like an image
   link in the body (``{attach}images/foo.png``, ``{static}/images/x.jpg``,
   an absolute URL, ...);
2. the first ``<img>`` in the rendered body (badge images, see
   ``THUMBNAIL_EXCLUDE_PATTERN``, are skipped);
3. the site-wide default ``THUMBNAIL_DEFAULT``.

Local images (and images embedded as ``data:`` URIs, as produced by
notebooks) are downscaled with Pillow into ``THUMBNAIL_PATH`` under the
output directory so that list pages stay light; if Pillow is not
installed the original image is used unchanged (a ``data:`` image is
then written out as a file as-is).  External URLs are used verbatim.
Identical sources yield one shared file (content-hashed names).

The result is stored on ``article.thumbnail`` as either a site-relative
path (to be prefixed with ``SITEURL`` by the template) or an absolute
URL; templates can tell them apart by the presence of ``//``.  The raw
metadata value, if any, is preserved as ``article.thumbnail_source``.

Settings:

- ``THUMBNAIL_DEFAULT``  : fallback image; a site-relative output URL
  (e.g. ``'images/logo.jpg'`` for a file under ``STATIC_PATHS``), an
  intrasite link (``'{static}/images/logo.jpg'``) or an absolute URL.
  ``None`` (default) means articles without an image get no thumbnail.
- ``THUMBNAIL_SIZE``     : bounding box for generated thumbnails, as
  (width, height); default ``(320, 240)``.  Images are only ever
  downscaled, keeping their aspect ratio.
- ``THUMBNAIL_PATH``     : output subdirectory for generated files
  (default ``'thumbnails'``).
- ``THUMBNAIL_QUALITY``  : JPEG quality for generated thumbnails
  (default 85).
- ``THUMBNAIL_EXCLUDE_PATTERN`` : regex; body images whose ``src``
  matches are never used (default ``r'shields\\.io'``, i.e. badges).
"""

import base64
import hashlib
import io
import os
import posixpath
import re
from urllib.parse import unquote, urlparse

from pelican import signals
from pelican.generators import ArticlesGenerator

import logging
logger = logging.getLogger(__name__)

try:
    from PIL import Image, ImageOps
except ImportError:  # pragma: no cover - optional dependency
    Image = None

_IMG_RE = re.compile(r'<img\b[^>]*?\bsrc\s*=\s*(["\'])(.*?)\1', re.I | re.S)
_INTRASITE_RE = re.compile(r'^\{(?P<what>static|attach|filename)\}(?P<value>.+)$')
_DATA_URI_RE = re.compile(r'^data:(?P<mime>image/[\w.+-]+)(?P<params>(;[\w.+-]+=[^;,]*)*)(?P<b64>;base64)?,(?P<data>.*)$', re.S)

_MIME_EXT = {'image/jpeg': '.jpg', 'image/png': '.png', 'image/gif': '.gif',
             'image/webp': '.webp', 'image/svg+xml': '.svg', 'image/bmp': '.bmp'}


class _Source:
    """An image candidate: where it can be read from and how it is addressed."""

    def __init__(self, url=None, path=None, data=None, ext=None):
        self.url = url        # site-relative output URL, or absolute URL
        self.path = path      # source file on disk (local images)
        self.data = data      # raw bytes (data: URIs)
        self.ext = ext        # file extension hint for data

    @property
    def external(self):
        return bool(self.url) and '//' in self.url


def _is_absolute_url(value):
    return value.startswith('//') or urlparse(value).scheme in ('http', 'https')


# ---------------------------------------------------------------- phase 1

def _intrasite_source_path(article, value):
    """Content-relative source path for an intrasite link value, as Pelican
    computes it (see ``Content.get_static_links``)."""
    path = urlparse(value).path
    if path.startswith('/'):
        path = path[1:]
    else:
        path = article.get_relative_source_path(
            os.path.join(article.relative_dir, path))
    return path.replace('%20', ' ')


def register_static_links(generator):
    """Make sure files named only in ``thumbnail`` metadata are picked up by
    the StaticGenerator (which runs after this signal)."""
    for article in generator.articles:
        spec = getattr(article, 'thumbnail', None)
        if not isinstance(spec, str):
            continue
        m = _INTRASITE_RE.match(spec.strip())
        if m and m.group('what') in ('static', 'attach'):
            generator.context['static_links'].add(
                _intrasite_source_path(article, m.group('value')))


# ---------------------------------------------------------------- phase 2

def _lookup_intrasite(article, what, value, context):
    key = 'generated_content' if what == 'filename' else 'static_content'
    store = context.get(key) or {}
    for candidate in (value, unquote(value)):
        obj = store.get(_intrasite_source_path(article, candidate))
        if obj is not None:
            return obj
    return None


def _resolve(spec, article, context):
    """Turn an image reference into a ``_Source`` (or None)."""
    spec = spec.strip()
    if not spec:
        return None

    m = _DATA_URI_RE.match(spec)
    if m:
        if not m.group('b64'):
            return None
        try:
            data = base64.b64decode(m.group('data'), validate=False)
        except (ValueError, TypeError):
            return None
        return _Source(data=data, ext=_MIME_EXT.get(m.group('mime'), '.bin'))

    if spec.startswith('data:'):
        return None

    m = _INTRASITE_RE.match(spec)
    if m:
        what, value = m.group('what'), m.group('value')
        obj = _lookup_intrasite(article, what, value, context)
        if obj is None:
            logger.warning("thumbnails: unable to find '%s' referenced by %s",
                           spec, article.get_relative_source_path())
            return None
        if what == 'attach':
            obj.attach_to(article)
        return _Source(url=obj.url, path=getattr(obj, 'source_path', None))

    if _is_absolute_url(spec):
        return _Source(url=spec)

    # A plain site-relative output URL.  For files copied verbatim from
    # STATIC_PATHS the URL equals the content-relative source path, so the
    # source file can be found (and resized) without touching the static
    # object's url/save_as -- reading those would freeze its output
    # location and break later {attach} relocation by Pelican.
    url = spec.lstrip('/')
    obj = (context.get('static_content') or {}).get(url)
    return _Source(url=url, path=getattr(obj, 'source_path', None))


def _first_body_image(article, exclude_re):
    content = getattr(article, '_content', None) or ''
    for m in _IMG_RE.finditer(content):
        src = m.group(2).strip()
        if not src or (exclude_re and exclude_re.search(src)):
            continue
        return src
    return None


class _Maker:
    """Creates (downscaled) thumbnail files under the output directory."""

    def __init__(self, settings):
        self.size = tuple(settings.get('THUMBNAIL_SIZE', (320, 240)))
        self.subdir = settings.get('THUMBNAIL_PATH', 'thumbnails').strip('/')
        self.quality = int(settings.get('THUMBNAIL_QUALITY', 85))
        self.output_path = settings['OUTPUT_PATH']
        self.cache = {}
        self.warned = False

    def make(self, source):
        """Return the site-relative URL to use for ``source``."""
        if source.external:
            return source.url
        if source.path is not None:
            try:
                with open(source.path, 'rb') as f:
                    data = f.read()
            except OSError:
                return source.url
            ext = os.path.splitext(source.path)[1].lower() or '.bin'
        elif source.data is not None:
            data, ext = source.data, source.ext
        else:
            return source.url

        key = hashlib.sha1(data).hexdigest()
        if key in self.cache:
            return self.cache[key]

        result = self._resize(data, ext)
        if result is None:
            # Not resizable (SVG, unknown format, or no Pillow).
            if source.url is not None:
                url = source.url
            else:
                url = self._write(key, ext, data)
        else:
            out_data, out_ext = result
            url = self._write(key, out_ext, out_data)
        self.cache[key] = url
        return url

    def _resize(self, data, ext):
        if Image is None:
            if not self.warned:
                logger.warning('thumbnails: Pillow is not installed; images '
                               'are used at their original size.')
                self.warned = True
            return None
        try:
            im = Image.open(io.BytesIO(data))
            im.load()
        except Exception:
            return None
        try:
            im = ImageOps.exif_transpose(im)
        except Exception:
            pass
        if im.mode in ('RGBA', 'LA', 'PA') or (
                im.mode == 'P' and 'transparency' in im.info):
            im = im.convert('RGBA')
            im.thumbnail(self.size, Image.LANCZOS)
            # Keep PNG only when transparency is actually used; images
            # that merely carry an opaque alpha channel (typical for
            # plotting libraries) are far smaller as JPEG.
            if im.getchannel('A').getextrema()[0] < 255:
                buf = io.BytesIO()
                im.save(buf, 'PNG', optimize=True)
                return buf.getvalue(), '.png'
        im = im.convert('RGB')
        im.thumbnail(self.size, Image.LANCZOS)
        buf = io.BytesIO()
        im.save(buf, 'JPEG', quality=self.quality, optimize=True, progressive=True)
        return buf.getvalue(), '.jpg'

    def _write(self, key, ext, data):
        name = key[:16] + ext
        url = posixpath.join(self.subdir, name)
        target = os.path.join(self.output_path, *url.split('/'))
        os.makedirs(os.path.dirname(target), exist_ok=True)
        if not os.path.exists(target) or os.path.getsize(target) != len(data):
            with open(target, 'wb') as f:
                f.write(data)
        return url


def add_thumbnails(generators):
    articles = None
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            settings = generator.settings
            context = generator.context
            articles = generator.articles
            break
    if not articles:
        return

    default = settings.get('THUMBNAIL_DEFAULT')
    pattern = settings.get('THUMBNAIL_EXCLUDE_PATTERN', r'shields\.io')
    exclude_re = re.compile(pattern) if pattern else None
    maker = _Maker(settings)

    for article in articles:
        spec = getattr(article, 'thumbnail', None)
        article.thumbnail_source = spec if isinstance(spec, str) else None
        article.thumbnail = None

        candidates = []
        if article.thumbnail_source:
            candidates.append(article.thumbnail_source)
        body = _first_body_image(article, exclude_re)
        if body:
            candidates.append(body)
        if default:
            candidates.append(default)

        for spec in candidates:
            source = _resolve(spec, article, context)
            if source is None:
                continue
            url = maker.make(source)
            if url:
                article.thumbnail = url
                break


def register():
    signals.article_generator_finalized.connect(register_static_links)
    signals.all_generators_finalized.connect(add_thumbnails)
