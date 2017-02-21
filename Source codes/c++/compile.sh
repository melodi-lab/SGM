#!/bin/bash

name="OUT_find_para_sup";
ls -hlt pipeline_v6_ch2_fast_ch2_super.cpp
ls -hlt sub_sup.cpp
ls -hlt find_para.cpp
g++ ./find_para.cpp ./pipeline.cpp ./sub_sup.cpp -O3 -o $name.out;



