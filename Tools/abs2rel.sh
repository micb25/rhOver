#!/bin/bash

# alpha, beta, energy, abs. energy (a.u.), rel. energy (a.u.), rel. energy (1/cm)

if [ -f "$1" ]; then
	MINVAL=`awk 'NR==1 {min=$3}; NR>1 && $3<min && $3>0 {min=$3}; END{print min}' "$1"`
	awk '{ if ($0 ~ /^$/) { print "" } else { printf "   %8.2f%8.2f%18.10f%18.10f%18.10f\n", $1, $2, $3, $3 - '$MINVAL', ($3 - '$MINVAL') * 219474.63633664 } }' "$1" > "$1.rel"
else
	echo ""
	echo "This script converts absolute electrostatic energies as obtained by a rhOver scan"
	echo "in the improved electrostatic approach into relative energies."
	echo ""
	echo -e "\tusage:\n\t\t$0 filename.dat"
fi

