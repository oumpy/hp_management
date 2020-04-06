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
elif [ `git branch -a | grep "remotes/origin/$sourcebranch"` ]; then
    $makecommand clean "OUTPUTDIR=./$outputdir/$previewdir/$sourcebranch" &&\
    echo "\nSITEURL += \'/$previewdir/$sourcebranch\'\n" >> ./content/contentpublishconf.py &&\
    $makecommand publish "OUTPUTDIR=./$outputdir/$previewdir/$sourcebranch"
else # branch deleted
    rm -rf "./$outputdir/$previewdir/$sourcebranch" 
fi &&\
cd $outputdir &&\
git add . &&\
git commit -m "${comment:-Update}" &&\
git push origin $targetbranch
