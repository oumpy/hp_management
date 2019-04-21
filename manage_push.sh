today=$(date "+%Y.%m.%d")

git add *
git commit -m "${today} commitment"
git push -u origin master
