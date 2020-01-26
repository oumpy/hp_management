today=$(date "+%Y.%m.%d")

export GLOBIGNORE='output'
git add *
git commit -m "${today} commitment"
git push -u origin master
