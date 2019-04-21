#! bin/bash

# https://qiita.com/driller/items/49a990cbdfb51afed620


# git clone https://github.com/yuna06/chinese.git
# cd english
# pelican-quickstart
pelican content -o output -s pelicanconf.py
make html
# make publish
ghp-import output
# git push origin gh-pages
git push https://github.com/oumpy/oumpy.github.io.git gh-pages:master
