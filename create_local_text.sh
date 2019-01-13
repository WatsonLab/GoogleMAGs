#!/bin/bash

mkdir -p runs

rm runs/*

for f in `cat fileofaccessions.txt`; do
       touch runs/$f.txt
done

