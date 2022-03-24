#!/bin/bash

for i in `ls ../instances/exact/exact01*` ; do
    ..//codes/clusteredit -i ${i}
done