#! /bin/bash

MEASUREMENTS=(1 2 3 4 5 6 7 10 20)
NAMES=('gaussian_by_column' 'gaussian_by_row')
OPTION=3

make clean
make

for NAME in ${NAMES[@]}; do
    SIZE=1

    echo "size,lu_time,ss_time" >> $NAME.csv

    for SIZE in ${MEASUREMENTS[@]}; do
      echo -n "$((SIZE))00," >> $NAME.csv
        ./ep2 $OPTION < testsfiles/m$SIZE.dat | awk -F': ' '$0=$2' | paste -sd "," >> $NAME.csv
        SIZE=$(( SIZE + 1 ))
    done

    OPTION=$(( OPTION + 1 ))
    mv $NAME.csv time-measurements/
done

