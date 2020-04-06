#!/bin/sh
comment="$1"
branch=${2:-master}
targetbranch=${3:-master}
lockid=$(sh ./3rdtools/misc-tools/pexlock -w -1 oumpy_repository) || exit 1 # do ex-lock
sh ./tools/gotobranch.sh $branch
sh ./tools/pushsite.sh "$comment" $targetbranch
sh ./3rdtools/misc-tools/punlock "$lockid"  # release the lock
