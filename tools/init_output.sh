#!/bin/sh
if [ ! -e "./output/.git" ]; then
    repository_url=`python3 -c "import os,sys; sys.path.append(os.curdir+'/content'); \
                    from contentpublishconf import SITEREPOSITORY; \
                    print(SITEREPOSITORY)"`
    rm -rf ./output
    git clone "$repository_url" ./output
fi
