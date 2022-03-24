#!/bin/bash

for i in `ls ../instances/exact/exact00*` ; do
    ..//codes/clusteredit -i ${i}
done