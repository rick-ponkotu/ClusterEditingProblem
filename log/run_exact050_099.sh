#!/bin/bash

for j in {5..9}; do
for i in `ls ../instances/exact/exact0${j}*` ; do
 timeout 1h ..//codes/clusteredit -i ${i} >> $1_exact.log || echo "${i} TLE" >> $1_exact.log 
done
done

