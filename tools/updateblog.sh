#!/bin/sh
comment="$1"
branch=${2:-master}
targetbranch=${3:-master}
sh ./tools/gotobranch.sh $branch
sh ./tools/pushblog.sh "$comment" $targetbranch
