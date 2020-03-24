#!/bin/sh
comment="$1"
branch=${2:-"master"}
makecommand=`which gmake`
makecommand=${makecommand:-`which make`}
cd output &&\
git pull origin $branch &&\
git checkout -f $branch &&\
cd ../ &&\
$makecommand clean &&\
$makecommand publish &&\
cd output &&\
git commit -a -m "${comment:-Update}" &&\
git push origin $branch
