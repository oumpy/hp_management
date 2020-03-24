#!/bin/sh
here=$(basename `pwd`)
if [ "$pwd" != "content" ]; then
  echo "Run this script in the 'content' directory."
  exit
fi
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
  directory="./pages"
else
  directory="./articles/${schoolyear}sy/$category"
fi
filepath="$directory/$filename.md"

if [ -e "$filepath" ]; then
  echo "$0 error: $filepath already exists."
else
  mkdir -p "$directory"
  while read line
  do
    echo $(eval echo $line)
  done < $template > "$filepath"
  echo "$0: $filepath created."
fi
