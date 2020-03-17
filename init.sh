#!/bin/sh
git clone https://github.com/oumpy/hp_management.git
cd hp_management
git clone https://github.com/oumpy/oumpy.github.io.git ./output
git submodule update -i
cd themes && git submodule update voidy-bootstrap