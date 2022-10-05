#!/bin/bash
  
FILENAME="missing_ft.txt"
in="/home/sareh/93/no_ovlp/hyphy_clean/"
out="/home/sareh/93/no_ovlp/missing_ft/"

LINES=$(cat $FILENAME)

for LINE in $LINES
do
    echo "$LINE"
    echo $out"$LINE"

    fasttree -gtr -nt -gamma $in"$LINE" > $out"$LINE"

done
