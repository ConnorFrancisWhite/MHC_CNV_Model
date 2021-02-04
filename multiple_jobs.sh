#!/bin/bash
for NUMBERS in {1..1}; do
	echo $NUMBERS
	./Poly ${NUMBERS} 5000 1000 10 11 3 0.0001 0.000001 0.0 0 0 1 0.0 0.01 0.1 100 1

done
