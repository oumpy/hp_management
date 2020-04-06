#!/bin/sh
comment="$1"
branch=${2:-master}
targetbranch=${3:-master}
lockid=$(pexlock -w -1 oumpy_repository) || exit 1 # do ex-lock
sh ./tools/gotobranch.sh $branch
sh ./tools/pushsite.sh "$comment" $targetbranch
punlock "$lockid"  # release the lock
