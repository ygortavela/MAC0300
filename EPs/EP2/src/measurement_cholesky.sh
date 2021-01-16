#! /bin/bash

MEASUREMENTS=(1 2 3 4 5 6 7 10 20)
NAMES=('cholesky_by_column' 'cholesky_by_row')
OPTION=1

make clean
make

for NAME in ${NAMES[@]}; do
    SIZE=1

    echo "size,chol_time,forw_time,back_time" >> $NAME.csv

    for SIZE in ${MEASUREMENTS[@]}; do
      echo -n "$((SIZE))00," >> $NAME.csv
        ./ep2 $OPTION < testsfiles/a$SIZE.dat | awk -F': ' '$0=$2' | paste -sd "," >> $NAME.csv
        SIZE=$(( SIZE + 1 ))
    done

    OPTION=$(( OPTION + 1 ))
    mv $NAME.csv time-measurements/
done

