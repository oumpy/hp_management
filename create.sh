#!/bin/sh
template="./template.md"
filename=$1
category=$2
tags=$3
date=$(date "+%Y.%m.%d")
year=$(date "+%Y")
month=$(date "+%m")
if [ $month -lt 4 ]; then
  schoolyear=$(($year-1))
else
  schoolyear=$year
fi

if [ "$category" = "" ] || [ "$category" = "page" ] || [ "$category" = "pages" ]; then
  directory="content/pages"
else
  directory="content/articles/${schoolyear}sy/$category"
  mkdir -p $directory
fi
while read line
do
  echo $(eval echo $line)
done < $template > "$directory/$filename.md"
