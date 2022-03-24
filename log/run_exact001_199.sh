#!/bin/bash

for k in {0..1}; do
for j in {0..9}; do
for i in `ls ../instances/exact/exact${k}${j}*` ; do
    ..//codes/clusteredit -i ${i}
done
done
done