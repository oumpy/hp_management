# AutoSummary Plugin for Pelican
## What is it ?
This plugin generates summary for articles & pages automatically,
by taking the beginning part of the contents.
Tags in the contents are kept, deleted or converted, belong to customizable settings.

## Setting Parameters & Default Values
- `AUTOSUMMARY_MAX_LENGTH = 140` : The maximum length (characters) of the generated summary.
- `AUTOSUMMARY_MIN_LENGTH = 1` : If the summary given in the metadata is shorter than this length, this plugin is applied.
- `AUTOSUMMARY_KEEP_HTMLTAGS = ('a', 'font', 's', 'strong', 'em', 'u', 'b')` : List of html tags which are preserved as it was.
- `AUTOSUMMARY_REPLACE_HTMLTAGS = {'h1':'strong', 'h2':'strong',..., 'h9':'strong'}` : Dictionary which indicates replacements of tags.

## Acqknowledgement & License
This plugin is forked from [Pelican Summary Plugin](https://github.com/getpelican/pelican-plugins/tree/master/summary).
Along with the upstream, this plugin is distributed under the GNU Affero General Public License v3.0 or later.
