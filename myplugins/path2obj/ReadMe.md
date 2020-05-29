# path2obj Plugin for Pelican
## What is it ?
This plugin provides fileters `source2obj` & `url2obj`,
which converts sourcefile paths or generated file paths
to the corresponding objects such as pages, article, tags, etc.

## How to Use
just set
```
PLUGINS += ['path2obj']
```
There are no setting parameters for this plugin.

In your theme, you can use it as
```
{% for source in AMAZING_SOURCE_LIST %}
    <p>{{ (source|source2obj).title }}</p>
{% endfor %}
```
where `AMAZING_SOURCE_LIST` is defined in `pelicanconf.py` as
```
PATH_LIST = ['article/amazing_announce.md', 'tag:amazing',...]
```

## License
This plugin is distributed under the GNU Affero General Public License version 3.0 (AGPLv3) or later.

## Author
This plugin was written by AtamaokaC <realatamaokac@gmail.com> & Osaka University Medical School Python Association (OUMPY).
