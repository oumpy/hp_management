#!/bin/sh
comment="$1"
branch=${2:-master}
targetbranch=${3:-master}
sh gotobranch.sh $branch
sh blog_push.sh "$comment" $targetbranch
