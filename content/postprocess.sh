#!/bin/sh
outputdir=${1:-output} 
cd $outputdir/articles/
ln -sfn ../blog/* ./
cd ../../
