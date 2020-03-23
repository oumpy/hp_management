#!/bin/sh
comment="$1"
branch=${2:-"master"}
git pull origin $branch
git checkout -f $branch
sh init.sh -c
sh blog_push $comment
