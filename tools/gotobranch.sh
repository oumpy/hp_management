#!/bin/sh
branch=${1:-"master"}
git fetch
if [ ! `git branch -a | sed 's/^[ \t]*//' | grep "^remotes/origin/$branch$"` ]; then
    git branch -D $branch
    echo "No branch $branch."
    exit
fi
rm -rf content/*
rm -rf theme/*
git checkout -f $branch
git pull origin $branch
sh ./tools/init.sh
