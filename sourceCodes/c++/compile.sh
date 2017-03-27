#!/bin/bash

name="OUT_SGM_sup";
ls -hlt pipeline.cpp
ls -hlt sub_sup.cpp
ls -hlt SGM.cpp
g++ ./SGM.cpp ./pipeline.cpp ./sub_sup.cpp -O3 -o $name.out;



