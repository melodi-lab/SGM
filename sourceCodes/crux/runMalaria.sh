#!/bin/bash

function crux_index {
rm -rf malaria-index
	$CRUX tide-index --output-dir $OUTPUT --clip-nterm-methionine T --enzyme lys-c --missed-cleavages $1 --peptide-list T \
	--mods-spec C+57.0214,K+229.16293 --nterm-peptide-mods-spec X+229.16293 \
	--max-mods 0 \
	--custom-enzyme '[K]|[X]' \
	    $TARGETS malaria-index
}
#CRUX="/n/trombone/s1/wrbai/codes/run_crux/crux-2.1.Linux.x86_64/bin/crux"
CRUX="/homes/wrbai/crux/crux-2.1.Linux.x86_64/bin/crux"
TARGETS="/n/trombone/s1/wrbai/codes/msgfplus/jeffH/plasmo_Pfalciparum3D7_NCBI.fasta"
#TARGETS="/n/trombone/s1/wrbai/codes/msgfplus/jeffH/test.fasta"
#P_TARGETS="/n/trombone/s1/wrbai/codes/msgfplus/jeffH/plasmo_Pfalciparum3D7_NCBI_dig.fasta"
OUTPUT="malaria-output"
#for i in 0 1 2 3 4 5 6 7 8 
#do
crux_index 8 
#done

