#!/bin/bash



MSGFPlusDir=$PWD


MSGFPlus="java -Xmx3500M -jar $MSGFPlusDir/MSGFPlus.jar"
MzIDToTSV="java -cp $MSGFPlusDir/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv"

OUTTARGETMZID="$outDir/human-ch$CHARGE-targets.mzid"
OUTTARGETTSV="$outDir/human-ch$CHARGE-targets.tsv"


$MzIDToTSV -i OUTTARGETMZID -o OUTTARGETTSV
	
