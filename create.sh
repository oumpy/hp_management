#!/bin/sh
template="./template.md"
filename=$1
category=$2
tags=$3
date=$(date "+%Y.%m.%d")
year=$(date "+%Y")

if [ "$category" = "" ] || [ "$category" = "page" ] || [ "$category" = "pages" ]; then
  directory="content/pages"
else
  directory="content/$year/$category"
  mkdir -p $directory
fi
while read line
do
  echo $(eval echo $line)
done < $template > "$directory/$filename.md"
