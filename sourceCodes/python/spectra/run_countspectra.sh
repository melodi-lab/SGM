#!/usr/bin/env bash
#
# Usage: run_countspectra.sh [ <directory> ]
#
# Prints out the number of spectra in each MS2 file in the specified
# directory. If no directory is given, default to the project data directory.

trap "{ echo \"$0: Script cancelled...\"; exit -2; }" SIGINT SIGTERM

DATADATA=$1

if [ ! -d "$DATADIR" ]; then
  DATADIR=/g/ssli/projects/msms/data/spectra
fi

for ms2file in $(ls $DATADIR/*.ms2*); do
  echo "Reading $ms2file"
  python countspectra.py -z -i $ms2file
done
