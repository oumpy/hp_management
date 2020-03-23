#!/bin/sh

# https://qiita.com/driller/items/49a990cbdfb51afed620


# git clone https://github.com/yuna06/chinese.git
# cd english
# pelican-quickstart
cd output
git pull origin master
cd -
# pelican content -o output -s pelicanconf.py
make clean
make publish
#ghp-import output
# git push origin gh-pages
# git push https://github.com/oumpy/oumpy.github.io.git master
cd output && git commit -a -m "Update" && git push origin master
