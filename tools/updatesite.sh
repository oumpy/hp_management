#!/bin/sh
comment="$1"
branch=${2:-master}
targetbranch=${3:-master}
sh ./tools/gotobranch.sh $branch
if [ "$branch"=="master" ]; then
    sh ./tools/setcgi.sh
fi
sh ./tools/pushsite.sh "$comment" $targetbranch
