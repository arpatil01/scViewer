#!/bin/bash

set -e
set -u

# Vendors knncolle and all the related libraries.

if [ ! -e source-knncolle ]
then 
    git clone https://github.com/LTLA/knncolle source-knncolle
else 
    cd source-knncolle
    git pull
    cd -
fi

cd source-knncolle
git checkout c5a1776ebf10641d9bf8715e89cb2d965b06e899
rm -rf ../knncolle
cp -r include/knncolle/ ../knncolle
git checkout master

cmake -S . -B build
for lib in aarand annoy kmeans powerit
do
    rm -rf ../${lib}
    cp -r build/_deps/${lib}-src/include/${lib} ../${lib}
done

rm -rf ../hnswlib
cp -r build/_deps/hnswlib-src/hnswlib ../hnswlib
cd -
