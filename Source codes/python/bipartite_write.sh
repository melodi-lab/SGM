#!/bin/bash

GENTESTDATA="write_msms_bipartite.py"
SHARDDATA="shard_spectra.py"
CURR=$(pwd)

TOL=3.0
# ORGANISM="worm-nokeil"
ORG="yeast"
PART=2
ORGANISM="$ORG-nokeil"
NORMALIZE=fastSequestTransform
# NORMALIZE=sequestProcess
MASS=2000
CHARGE=3
ENCODE="$CURR/$ORG-02-ch$CHARGE"
MS2="/s1/wrbai/codes/ms2file/${ORG}_parts"

if [ ! -d $ENCODE ]
then
    mkdir $ENCODE
else
    rm -rfv $ENCODE
    mkdir $ENCODE
fi

echo "charge=$CHARGE"

#MS2FILE="test.ms2"
MS2FILE="$MS2/$ORG-0$PART-ch$CHARGE.ms2"
echo "Converting file $MS2FILE"
python $SHARDDATA \
    --spectra $MS2FILE \
    --num_spectra $(grep -c '^S' $MS2FILE) \
    --targets $ORGANISM.fasta.peptides \
    --decoys $ORGANISM-decoy.fasta.peptides \
    --tol $TOL \
    --shards 1 \
    --output_dir $ENCODE \
    --charges $CHARGE

for pickle in $ENCODE/shard-*.pickle
do
    python $GENTESTDATA --output $ENCODE --pickle $pickle \
	--max_mass 2000 --targets_per_spectrum 100000 \
	--normalize $NORMALIZE  \
	--bin_width 1.0005079 --bin_offset 0.68 --charge $CHARGE \
	--spec_file "$ORG-0$PART-ch$CHARGE" --bipartite_file "msms-$NORMALIZE-bw1.0005-bo0.68-tol3Da-$ORG-0$PART-ch$CHARGE.txt" \
	--xcorr_all_psms_file "xcorrAllPSMs-$ORG-0$PART-ch$CHARGE.txt" \
	--xcorr_ident_file "xcorrIdent-$ORG-0$PART-ch$CHARGE.txt" \
	--do_not_write_all_psms
done
