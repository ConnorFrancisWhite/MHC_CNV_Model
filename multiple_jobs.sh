#!/bin/bash
for NUMBERS in {52..52}; do
	echo $NUMBERS
	./Poly ${NUMBERS} 5000 50000 10 21 3 0.0001 0.000001 0.0 1 0 0 0.0 10000

done
