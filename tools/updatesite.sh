#!/bin/sh
comment="$1"
sourcebranch=${2:-master}
targetbranch=${3:-master}
lockid=$(sh ./3rdtools/misc-tools/pexlock -w -1 oumpy_repository) || exit 1 # do ex-lock
sh ./tools/gotobranch.sh $sourcebranch
if [ "$sourcebranch" = "master" ]; then
    sh ./tools/setcgi.sh
fi
sh ./tools/pushsite.sh "$comment" $sourcebranch $targetbranch
sh ./3rdtools/misc-tools/punlock "$lockid"  # release the lock
