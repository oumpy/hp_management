# Subsections Plugin for Pelican
## What is it ?
This plugin adds `subsections` attribute to each article/page,
which contain tree-structure information of the sections hierarchy.

## How to Use
Just add `subsections` to `PLUGINS` in your `pelicanconf.py`.

The `subsections` is a `list` of `Sections` object, whose attributes are:
| attribute | type | description | example |
|---|---|---|---|
| `title` | `str` | title of the section | `Introduction` |
| `url` | `str` | url with the id of the section. | `articles/great.html#pelican-subsections-7` |
| `subsections` | `list` of `Sections` | child-sections (grandchildren not included) | |

For example, if the content of an article or a page `obj` is like
```
<h2>Introduction</h2>
<h4>Surprising</h4>
<h5>Awesome</h5>
<h3>Amazing</h3>
<h2>Discussion</h2>
<h1>Conclusion</h1>
...
```
then this plugin creates `obj.subsections` as
```
obj.subsection = [
    Section('Introduction', 'articles/great.html#pelican-subsections-0', [
        Section('Surprising', 'articles/great.html#pelican-subsections-1', [
            Section('Awesome', 'articles//great.html#pelican-subsections-2', []),
        ]),
        Section('Amazing', 'articles/great.html#pelican-subsections-3', []),
    ]),
    Section('Discussion', 'articles/great.html#pelican-subsections-4', []),
    Section('Conclusion', 'articles/great.html#pelican-subsections-5', []),
    ...
]
```
where `<span id="pelican-subsections-xx">...</span>` are automatically inserted just inside  `<hx>...</hx>` in the content. 

## License
This plugin is distributed under the GNU Affero General Public License version 3.0 (AGPLv3) or later.

## Author
This plugin was written by AtamaokaC <realatamaokac@gmail.com> & Osaka University Medical School Python Association (OUMPY).

