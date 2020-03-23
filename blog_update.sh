#!/bin/sh
comment="$1"
branch=${2:-"master"}
targetbranch=${3:-"master"}
rm -rf content/*
rm -rf theme/*
git checkout -f $branch
git pull origin $branch
sh init.sh -c
sh blog_push.sh "$comment" $targetbranch
