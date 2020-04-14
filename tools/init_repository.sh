#!/bin/sh
git submodule update -i
cd plugins
git submodule update --init pelican-ipynb
cd ..

themename="voidy-bootstrap"
cd themes
git submodule update --init "$themename"
cd ..
mkdir -p "theme/$themename"
cp -an "themes/$themename"/* "theme/$themename/"
