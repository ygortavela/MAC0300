#! /bin/bash

MEASUREMENTS=11
NAMES=('matrix_matrix_jki' 'matrix_matrix_ikj')
OPTION=5

make clean
make time

for NAME in ${NAMES[@]}; do
    i=1
    SIZE=2

    while [ "$i" -le $MEASUREMENTS ]; do
        echo -n "$SIZE," >> $NAME.log
        cat matrix-test/matrix-$SIZE matrix-test/matrix-$SIZE | ./ep1 $OPTION | awk -F': ' '$0=$2' >> $NAME.log
        SIZE=$(( SIZE * 2 ))
        i=$(( i + 1 ))
    done

    OPTION=$(( OPTION + 1 ))
    mv $NAME.log time-results/
done

