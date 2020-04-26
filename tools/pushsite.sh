#!/bin/sh
comment="$1"
sourcebranch=${2:-master}
targetbranch=${3:-master}
makecommand=`which gmake`
makecommand=${makecommand:-`which make`}
previewdir="preview"
outputdir="output"
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
    $makecommand clean &&\
    mv $previewdir output/ &&\
    $makecommand publish
else
    git fetch
    if [ `git branch -a | sed 's/^[ \t]*//' | grep "^remotes/origin/$sourcebranch$"` ]; then
        git branch -D $sourcebranch
        $makecommand clean "OUTPUTDIR=./$outputdir/$previewdir/$sourcebranch" &&\
        {
            echo
            echo "SITEURL += '/$previewdir/$sourcebranch'"
            echo "if 'PREVIEW_SITENAME_APPEND' in globals():"
            echo "    SITENAME += PREVIEW_SITENAME_APPEND"
            echo "    SITETAG += PREVIEW_SITENAME_APPEND"
            echo
        }  >> ./content/contentconf.py &&\
        $makecommand html "OUTPUTDIR=./$outputdir/$previewdir/$sourcebranch"
    else # branch deleted
        rm -rf "./$outputdir/$previewdir/$sourcebranch" 
    fi
fi &&\
cd $outputdir &&\
git add . &&\
git commit -m "${comment:-Update}" &&\
git push origin $targetbranch
