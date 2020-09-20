#! /bin/bash

MEASUREMENTS=13
NAMES=('matrix_vector_by_row' 'matrix_vector_by_column')
OPTION=3

make clean
make time

for NAME in ${NAMES[@]}; do
    i=1
    SIZE=2

    while [ "$i" -le $MEASUREMENTS ]; do
        echo -n "$SIZE," >> $NAME.log
        cat matrix-test/matrix-$SIZE vector-test/vector-$SIZE | ./ep1 $OPTION | awk -F': ' '$0=$2' >> $NAME.log
        SIZE=$(( SIZE * 2 ))
        i=$(( i + 1 ))
    done

    OPTION=$(( OPTION + 1 ))
    mv $NAME.log time-results/
done

