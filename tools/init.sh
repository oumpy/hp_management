#!/bin/sh
if [ ! -e "./output/.git" ]; then
    repository_url=`python3 -c "import os,sys; sys.path.append(os.curdir+'/content'); \
                    from contentpublishconf import SITEREPOSITORY; \
                    print(SITEREPOSITORY)"`
    rm -rf ./output
    git clone "$repository_url" ./output
fi

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
