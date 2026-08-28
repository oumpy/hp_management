#!/bin/sh
comment="$1"
sourcebranch=${2:-master}
targetbranch=${3:-master}
previewdir="preview"
outputdir="output"
# Empty an output directory while keeping dotfiles (in particular .git).
clean_output () {
    [ ! -d "$1" ] || rm -rf "$1"/*
}
cd $outputdir &&\
git pull origin $targetbranch &&\
git checkout -f $targetbranch &&\
cd ../ &&\
if [ "$sourcebranch" = "master" ]; then
    if [ -d "./$outputdir/$previewdir" ]; then
        rm -rf $previewdir && mv "./$outputdir/$previewdir" ./
    else
        mkdir $previewdir
    fi &&\
    clean_output "$outputdir" &&\
    mv $previewdir output/ &&\
    pelican -s publishconf.py
else
    git fetch
    if [ `git branch -a | sed 's/^[ \t]*//' | grep "^remotes/origin/$sourcebranch$"` ]; then
        git branch -D $sourcebranch
        clean_output "./$outputdir/$previewdir/$sourcebranch" &&\
        {
            echo
            echo "SITEURL += '/$previewdir/$sourcebranch'"
            echo "if 'PREVIEW_SITENAME_APPEND' in globals():"
            echo "    SITENAME += PREVIEW_SITENAME_APPEND"
            echo "    SITETAG += PREVIEW_SITENAME_APPEND"
            echo
        }  >> ./content/contentconf.py &&\
        pelican -o "./$outputdir/$previewdir/$sourcebranch"
    else # branch deleted
        rm -rf "./$outputdir/$previewdir/$sourcebranch" 
    fi
fi &&\
cd $outputdir &&\
git add . &&\
git commit -m "${comment:-Update}" &&\
git push origin $targetbranch
