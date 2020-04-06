#!/bin/sh
outputdir=${1:-output} 
cd $outputdir/articles/
ln -s ../blog/* ./
cd ../../
