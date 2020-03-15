#!/bin/sh
template="./template.md"
filename=$1
category=$2
tags=$3
if [ "$category" = "" ]; then
  category="pages"
fi
mkdir -p content/$category
date=$(date "+%Y.%m.%d")
while read line
do
  echo $(eval echo $line)
done < $template > content/$category/$filename.md
