# RenameCategory Plugin for Pelican
## What is it ?
Using this plugin, you can change or add names of categories.
The slugs are not changed, so the URLs for the category pages are kept as they are.

## Customizable Parameters
`RENAMECATEGORY_ATTRIBUTES = ('name', 'name_long')` : attributes for which new names are assigned
`RENAMECATEGORY_ALTERNATIVES = dict()` : Dictionary of `(category, list of new names)`

## How to Use
For example, if you set in your `pelicanconf.py` as
```
RENAMECATEGORY_ATTRIBUTES = ('name', 'name_long')
RENAMECATEGORY_ALTERNATIVES = { 'beginners' : ('ABC', 'AtCoder Beginner Contest').
                                'grand' : ('AGC', 'AtCoderGrand Contest') }
```
a statement `{{ catetory }}` or `{{ article.category }}` in your theme prints `'ABC'` for category `'beginners'`.
This can easily change the appearance without changing the theme at all.
Furthermore, by putting `{{ category.name_long }}` in your theme, you can obtain `'AtCoder Beginner Contest'`.

## License
This plugin is distributed under the GNU Affero General Public License version 3.0 (AGPLv3) or later.

## Author
This plugin was written by AtamaokaC <realatamaokac@gmail.com> & Osaka University Medical School Python Association (OUMPY).
