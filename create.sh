#! bin/bash
touch content/$1.md
today=$(date "+%Y.%m.%d")
echo "Title:
Date: ${today}
Category: $2
Tags: $2
Slug: $1
Author: 
Summary:
" > content/$1.md
open content/$1.md
cd content
mv $1.md $2/$1.md
cd ../
