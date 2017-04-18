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
	-ti 0,0 \
	-mod jeffH/mods.txt \
	-minLength 7 \
	-thread 1 \
	-ntt 2 \
	-ignoreMetCleavage 1 
	#> $3.log.txt
	# -mod msgf-mods.txt \

    $MzIDToTSV -i $3 \
        -o $5 
}

MSGFPlusDir="/homes/wrbai/co/"


MSGFPlus="java -Xmx3500M -jar $MSGFPlusDir/MSGFPlus.jar"
MzIDToTSV="java -cp $MSGFPlusDir/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv"
TARGETDB="/n/trombone/s1/wrbai/codes/run_crux/Linfeng-output-2/tide-index.peptides.target.fasta"
DECOYDB="/n/trombone/s1/wrbai/codes/run_crux/Linfeng-output-2/tide-index.peptides.decoy.fasta"
#MS2FILE="/n/trombone/s1/wrbai/codes/ms2file/Linfeng/Linfeng_010211_HapMap32_1.ms2"

for (( i=1; i<=20; i++))
do
	id=$((i*5-5+$1))
	IID=$(($1)) 
	#MS2FILE="/n/trombone/s0/wrbai/codes/ms2file/Linfeng/Linfeng_080510_HapMap9_$id.ms2"
	for(( c=2; c<=5; c++))
	do
	
	    CHARGE=$c
	    outDir="linfeng-$id/linfeng-ch$CHARGE"
	    #outDir="linfeng_112710/linfeng-ch$CHARGE"
		#outDir="linfeng_9_$IID-1/linfeng-ch$CHARGE"
		echo $outDir
	    if [[ ! -d $outDir ]]
	    then
		mkdir -p $outDir
	    else
		rm -rf $outDir
		mkdir -p $outDir
	    fi
	
	#    MS2="/n/trombone/s0/wrbai/codes/ms2file/Linfeng/Linfeng_080510_HapMap9_$IID-ch$c.mgf"
	    MS2="/n/trombone/s0/wrbai/codes/ms2file/Linfeng/all/Linfeng-$id-ch$c.mgf"
	#    MS2="/s1/wrbai/program/Linfeng_010211_HapMap32_1-ch3.mgf"
	#	MS2="../ms2file/malaria/malaria_TMT_11-ch3.mgf"
		echo $MS2
	 #   MS2="/n/trombone/s1/wrbai/codes/ms2file/highRes/LinfengPrunedZeroIntensities_112710_HapMap27_1-ch$c.mgf"
	# start with targets
	# start with targets
	    OUTTARGETMZID="$outDir/linfeng-ch$CHARGE-targets.mzid"
	    OUTTARGETTSV="$outDir/linfeng-ch$CHARGE-targets.tsv"
	    runMsgf $MS2 $TARGETDB $OUTTARGETMZID $CHARGE $OUTTARGETTSV
	
	# next decoys
	    OUTDECOYMZID="$outDir/linfeng-ch$CHARGE-decoys.mzid"
	    OUTDECOYTSV="$outDir/linfeng-ch$CHARGE-decoys.tsv"
	    runMsgf $MS2 $DECOYDB $OUTDECOYMZID $CHARGE $OUTDECOYTSV
	
	# parse results
	#    IDENT="$outDir/linfeng-ch$CHARGE-ident.txt"
	 #   echo "$OUTTARGETTSV $OUTDECOYTSV"
	  #  parse $OUTTARGETTSV $OUTDECOYTSV $MS2FILE $CHARGE $IDENT
	done
done
