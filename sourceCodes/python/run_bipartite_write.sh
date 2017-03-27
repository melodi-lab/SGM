#!/bin/bash

CHARGE="3"
ORG="plasm-1"
MS2FILE="data/$ORG-ch$CHARGE.ms2"
TARGET="data/tide-index.peptides.target.txt" 
DECOY="data/tide-index.peptides.decoy.txt"
./bipartite_write.sh $ORG $CHARGE $MS2FILE $TARGET $DECOY
