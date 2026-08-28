#!/bin/sh
# Initialize submodules required for site management.
# (The theme and all pelican plugins are now vendored in the repository;
#  only 3rdtools/misc-tools, used by legacy webhook deployment scripts,
#  remains as a submodule.)
git submodule update --init 3rdtools/misc-tools
