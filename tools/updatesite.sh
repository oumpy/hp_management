#!/bin/sh
comment="$1"
sourcebranch=${2:-master}
targetbranch=${3:-master}
sh ./tools/gotobranch.sh $sourcebranch
if [ "$sourcebranch" = "master" ]; then
    sh ./tools/setcgi.sh
fi
sh ./tools/pushsite.sh "$comment" $sourcebranch $targetbranch
