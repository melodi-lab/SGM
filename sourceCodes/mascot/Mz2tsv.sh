#!/bin/bash



MSGFPlusDir=$PWD


MzIDToTSV="java -cp $MSGFPlusDir/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv"

OUTTARGETMZID="test/F004917.mzid"
OUTTARGETTSV="test/F004917.tsv"


$MzIDToTSV -i $OUTTARGETMZID -o $OUTTARGETTSV
	
