#!/bin/sh
branch=${1:-"master"}
rm -rf content/*
rm -rf theme/*
git checkout -f $branch
git pull origin $branch
sh init.sh -c
