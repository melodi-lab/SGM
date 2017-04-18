#!/bin/bash

function parse {
    ./parseMsgfPlus.py --msgfPlusTargetFile $1 \
	--msgfPlusDecoyFile $2 \
	--spectra $3 \
	--charge $4 \
	--ident $5
}

function runMsgf {
    $MSGFPlus \
	-s $1 \
	-d $2 \
	-o $3 \
	-t 50ppm \
	-tda 0 \
	-n 1 \
	-minCharge $4 \
	-maxCharge $4 \
	-m 1 \
	-inst 1 \
	-e 1 \
	-ti 0,0 \
	-minLength 7 \
	-mod tmtModFile.txt \
	-thread 1 \
	-ntt 2 \
	-ignoreMetCleavage 1 
	#> $3.log.txt
	#-mod msgf-mods.txt \

    $MzIDToTSV -i $3 \
        -o $5 
}

MSGFPlusDir=$PWD


MSGFPlus="java -Xmx3500M -jar $MSGFPlusDir/MSGFPlus.jar"
MzIDToTSV="java -cp $MSGFPlusDir/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv"

TARGETDB="../data/tide-index.peptides.target.fasta"
DECOYDB="../data/tide-index.peptides.decoy.fasta"

	id=$1
c=$2	
	    CHARGE=$c
	    outDir="plasm-no-isotopy/plasm-$id/plasm-ch$CHARGE"
	echo $outDir
	    if [[ ! -d $outDir ]]
	    then
		mkdir -p $outDir
	    else
		rm -rf $outDir
		mkdir -p $outDir
	    fi
	
		MS2="../data/plasm-$id-ch$c.mgf"
		echo $MS2
	# start with targets
	    OUTTARGETMZID="$outDir/plasm-ch$CHARGE-targets.mzid"
	    OUTTARGETTSV="$outDir/plasm-ch$CHARGE-targets.tsv"
	    runMsgf $MS2 $TARGETDB $OUTTARGETMZID $CHARGE $OUTTARGETTSV
	
	# next decoys
	    OUTDECOYMZID="$outDir/plasm-ch$CHARGE-decoys.mzid"
	    OUTDECOYTSV="$outDir/plasm-ch$CHARGE-decoys.tsv"
	    runMsgf $MS2 $DECOYDB $OUTDECOYMZID $CHARGE $OUTDECOYTSV
	
