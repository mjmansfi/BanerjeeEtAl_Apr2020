#!/bin/bash

PROTS=("ACE2" "DPP4")
ANIMALS=("mammals" "others")

for PROT in ${PROTS[@]};
do
	for ANIMAL in ${ANIMALS[@]};
	do
		# if output already exists, delete it
		if [[ -f "$PROT.ordered.shared.motif_$ANIMAL.tsv" ]];
		then
			rm $PROT.ordered.shared.motif_$ANIMAL.tsv
		fi

		grep ">" $PROT.ordered.shared.motif_$ANIMAL.mfa | cut -d ">" -f 2 | cut -d '/' -f 1 | while read ID;
			do
				SEQ=$(grep -A 1 "$ID" $PROT.ordered.shared.motif_others.mfa | tail -n 1);
				#printf "$ID\t$SEQ\t$PROT\t$ANIMAL\n"
				printf "$ID\t$SEQ\t$PROT\t$ANIMAL\n" >> "$PROT.ordered.shared.motif_$ANIMAL.tsv"
			done
	done
done
#grep ">" ACE2.ordered.shared.motif_others.mfa
# | cut -d ">" -f 2 |
# cut -d '/' -f 1 |
# while read ID;
# do
# SEQ=$(grep -A 1 "$ID" ACE2.ordered.shared.motif_others.mfa | tail -n 1);
# printf "$ID\t$SEQ\tNon-mammal\tACE2\n";
# done > ACE2.ordered.shared.motif_others.^
