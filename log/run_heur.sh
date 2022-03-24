#!/bin/bash

for j in {0..2}; do
for i in `ls ../instances/exact/exact0${j}*` ; do
 timeout 1m ..//codes/clusteredit -i ${i} >> $1_heur.log || echo "${i} TLE" >> $1_heur.log
done
done