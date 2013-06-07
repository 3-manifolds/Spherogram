#!/bin/bash
svn checkout http://planarity.googlecode.com/svn/trunk/ planarity-read-only
cd planarity-read-only
patch -p0 < ../planarity.patch
cd c
if [ `uname` = "Darwin" ]; then
  gcc -arch i386 -arch x86_64 -c *.c nauty/*.c
  gcc -arch i386 -arch x86_64 -o planarity *.o
else
  gcc -c -fPIC *.c nauty/*.c
  gcc -o planarity *.o
fi