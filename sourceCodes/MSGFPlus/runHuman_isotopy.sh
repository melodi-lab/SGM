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
	-t 10ppm \
	-tda 0 \
	-n 1 \
	-minCharge $4 \
	-maxCharge $4 \
	-m 1 \
	-inst 1 \
	-e 1 \
	-ti 0,1 \
	-minLength 7 \
	-mod mods.txt \
	-thread 1 \
	-ignoreMetCleavage 1 
	#> $3.log.txt
	#-mod msgf-mods.txt \

    $MzIDToTSV -i $3 \
        -o $5 
}

MSGFPlusDir=$PWD


MSGFPlus="java -Xmx3500M -jar $MSGFPlusDir/MSGFPlus.jar"
MzIDToTSV="java -cp $MSGFPlusDir/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv"

TARGETDB="database/ipi.HUMAN.v3.74-nokeil-target.fasta"
DECOYDB="database/ipi.HUMAN.v3.74-nokeil-decoy.fasta"

	id=$1
c=$2	
	    CHARGE=$c
	    outDir="human-with-isotopy/human-$id/human-ch$CHARGE"
	echo $outDir
		MS2="/home/ubuntu/SGM/sourceCodes/data/msgfdata/human/Linfeng-$id-ch$c.mgf"
if [ ! -f $MS2 ]; then
	exit  
fi
	    if [[ ! -d $outDir ]]
	    then
		mkdir -p $outDir
	    else
		rm -rf $outDir
		mkdir -p $outDir
	    fi
	
		echo $MS2
	# start with targets
	    OUTTARGETMZID="$outDir/human-ch$CHARGE-targets.mzid"
	    OUTTARGETTSV="$outDir/human-ch$CHARGE-targets.tsv"
	    runMsgf $MS2 $TARGETDB $OUTTARGETMZID $CHARGE $OUTTARGETTSV
	
	# next decoys
	    OUTDECOYMZID="$outDir/human-ch$CHARGE-decoys.mzid"
	    OUTDECOYTSV="$outDir/human-ch$CHARGE-decoys.tsv"
	    runMsgf $MS2 $DECOYDB $OUTDECOYMZID $CHARGE $OUTDECOYTSV
	
