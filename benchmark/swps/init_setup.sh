#!/bin/sh
swps_dir=$(dirname `realpath $0`)
cd $swps_dir

wget http://lab.dessimoz.org/swps3/files/swps3-src-current.tar.bz2
tar xjvf swps3-src-current.tar.bz2

rm swps3-src-current.tar.bz2
cd swps3-*

wget http://lab.dessimoz.org/swps3/files/matrices.tar.bz2
tar xjvf matrices.tar.bz2

wget http://lab.dessimoz.org/swps3/files/test.tar.bz2
tar xjvf test.tar.bz2 
