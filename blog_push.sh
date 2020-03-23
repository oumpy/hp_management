#!/bin/sh
cd output && git pull origin master && cd -
make clean
make publish
cd output && git commit -a -m "Update" && git push origin master
