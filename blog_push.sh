#!/bin/sh
comment="$1"
cd output && git pull origin master && cd -
make clean
make publish
cd output && git commit -a -m ${comment:-"Update"} && git push origin master
