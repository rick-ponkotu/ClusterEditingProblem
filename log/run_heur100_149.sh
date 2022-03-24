#!/bin/bash

for j in {0..4}; do
for i in `ls ../instances/heur/heur1${j}*` ; do
 timeout 1h ..//codes/clusteredit -i ${i} >> $1_heur.log || echo "${i} TLE" >> $1_heur.log 
done
done

