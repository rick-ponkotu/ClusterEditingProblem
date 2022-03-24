#!/bin/bash

for i in `ls ../instances/exact/*` ; do
 j=$(head -n 1 ${i})
 echo "${i} ${j}" | cut -c 20-30,37-
done
for i in `ls ../instances/heur/*` ; do
 j=$(head -n 1 ${i})
 echo "${i} ${j}"  | cut -c 19-28,35- 
done
