#!/bin/bash
#
# How the source was downloaded.  

svn checkout http://planarity.googlecode.com/svn/trunk/ planarity-read-only
cd planarity-read-only
patch -p0 < ../planarity.patch
cd ../
mv planarity-read-only/c .
rm -r  planarity-read-only
rm -r  c/.cproject c/.project c/.settings/ c/*.orig c/samples