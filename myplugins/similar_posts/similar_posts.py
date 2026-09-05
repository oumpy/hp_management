# -*- coding: utf-8 -*-
"""
SimilarPosts Plugin for Pelican (in-house, content-based)
---------------------------------------------------------

Selects related articles by *content* similarity instead of the
shared-tag counting of the former pelican-related-posts plugin.  Tag
counting has little discriminative power when most articles carry only
a few tags drawn from a small set of large buckets.

The plugin is named ``similar_posts`` to avoid confusion with the
original; the attribute it sets (``article.related_posts``) and the
``RELATED_*`` setting names are kept for template compatibility.

Method (pure Python, no new dependencies):

- Each article's text (title, weighted higher, plus HTML-stripped body)
  is represented as a TF-IDF vector over character n-grams
  (RELATED_NGRAM_SIZES, default n = 2 and 3).  Character n-grams work
  well for Japanese without a morphological analyzer, and handle mixed
  Japanese/English/code text gracefully.
- Similarity is the cosine between vectors, plus a small bonus for
  shared tags weighted by tag rarity (IDF), so tags still contribute
  where they are informative.
- The top RELATED_POSTS_MAX articles (default 3) whose score reaches
  RELATED_MIN_SCORE are stored on ``article.related_posts`` -- the same
  attribute the theme template already uses, so no template change is
  needed.  The list may therefore be shorter than the maximum, or empty
  for an article with no sufficiently similar peer; templates should
  treat the attribute as a plain, possibly empty list.
- Deterministic: ties are broken by date (newer first), then slug.

Settings:

- ``RELATED_POSTS_MAX``  : number of related articles (default 3).
- ``RELATED_MIN_SCORE``  : minimum similarity score for an article to
  be listed at all (default 0.0, i.e. no cutoff).  Raising it trades
  list completeness for precision: weakly related filler entries are
  dropped instead of being shown.
- ``RELATED_NGRAM_SIZES``: character n-gram lengths used for the
  content vectors (iterable of ints, default ``(2, 3)``).  Note that
  scores shift with this choice, so ``RELATED_MIN_SCORE`` may need
  retuning when it changes.
- ``RELATED_TAG_WEIGHT`` : weight of the tag-similarity bonus added to
  the content cosine (default 0.2).
- ``RELATED_TEXT_LIMIT`` : max characters of body text considered per
  article (default 20000; bounds the influence of huge notebooks).
"""

import html as html_module
import math
import re
import unicodedata
from collections import Counter

from pelican import signals
from pelican.generators import ArticlesGenerator

import logging
logger = logging.getLogger(__name__)

_DATA_URI_RE = re.compile(r'(src|href)="data:[^"]*"')
_SCRIPT_RE = re.compile(r'<(script|style)\b.*?</\1>', re.S | re.I)
_TAG_RE = re.compile(r'<[^>]+>')
_SPACE_RE = re.compile(r'\s+')


def _extract_text(article, limit):
    text = article._content or ''
    text = _DATA_URI_RE.sub(' ', text)
    text = _SCRIPT_RE.sub(' ', text)
    text = _TAG_RE.sub(' ', text)
    text = html_module.unescape(text)
    text = _SPACE_RE.sub(' ', text)
    text = unicodedata.normalize('NFKC', text).lower()[:limit]
    title = unicodedata.normalize('NFKC', article.title or '').lower()
    # Weight the title by repeating it: it is the strongest topical signal.
    return (title + ' ') * 5 + text


def _ngrams(text, sizes):
    counts = Counter()
    for n in sizes:
        for i in range(len(text) - n + 1):
            gram = text[i:i + n]
            if gram.isspace():
                continue
            counts[gram] += 1
    return counts


def _tfidf_vector(counts, df, n_docs):
    vec = {}
    for gram, tf in counts.items():
        # Grams appearing in only one document cannot contribute to any
        # cosine between different documents; dropping them shrinks the
        # vectors (and build time) considerably at negligible cost.
        if df[gram] < 2:
            continue
        idf = math.log((1 + n_docs) / (1 + df[gram])) + 1.0
        vec[gram] = (1.0 + math.log(tf)) * idf
    norm = math.sqrt(sum(w * w for w in vec.values()))
    if norm > 0:
        for gram in vec:
            vec[gram] /= norm
    return vec


def _cosine(v1, v2):
    if len(v2) < len(v1):
        v1, v2 = v2, v1
    return sum(w * v2[g] for g, w in v1.items() if g in v2)


def add_related_posts(generators):
    articles = None
    for generator in generators:
        if isinstance(generator, ArticlesGenerator):
            settings = generator.settings
            articles = generator.articles
            break
    if not articles:
        return

    max_count = settings.get('RELATED_POSTS_MAX', 3)
    min_score = settings.get('RELATED_MIN_SCORE', 0.0)
    ngram_sizes = tuple(settings.get('RELATED_NGRAM_SIZES', (2, 3)))
    tag_weight = settings.get('RELATED_TAG_WEIGHT', 0.2)
    text_limit = settings.get('RELATED_TEXT_LIMIT', 20000)
    n_docs = len(articles)

    # --- content vectors (TF-IDF over character n-grams)
    counts_list = [_ngrams(_extract_text(a, text_limit), ngram_sizes)
                   for a in articles]
    df = Counter()
    for counts in counts_list:
        df.update(counts.keys())
    vectors = [_tfidf_vector(c, df, n_docs) for c in counts_list]

    # --- tag IDF (rare tags say more than huge buckets)
    tags_list = [frozenset(str(t) for t in getattr(a, 'tags', []) or [])
                 for a in articles]
    tag_df = Counter(t for tags in tags_list for t in tags)
    tag_idf = {t: math.log(n_docs / c) for t, c in tag_df.items()}
    max_tag_idf = max(tag_idf.values()) if tag_idf else 1.0

    def tag_sim(i, j):
        shared = tags_list[i] & tags_list[j]
        if not shared or max_tag_idf == 0:
            return 0.0
        return min(1.0, sum(tag_idf[t] for t in shared) / max_tag_idf)

    # --- score all pairs, pick top-k per article
    for i, article in enumerate(articles):
        scored = []
        for j, other in enumerate(articles):
            if i == j:
                continue
            score = _cosine(vectors[i], vectors[j]) + tag_weight * tag_sim(i, j)
            if score < min_score:
                continue
            scored.append((-score, other.date, other.slug, other))
        scored.sort(key=lambda x: (x[0], -x[1].toordinal(), x[2]))
        article.related_posts = [s[3] for s in scored[:max_count]]


def register():
    signals.all_generators_finalized.connect(add_related_posts)
