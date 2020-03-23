#!/bin/sh
comment="$1"
branch=${2:-"master"}
cd output &&\
git pull origin $branch &&\
git checkout -f $branch &&\
cd ../ &&\
make clean &&\
make publish &&\
cd output &&\
git commit -a -m ${comment:-"Update"} &&\
git push origin $branch
