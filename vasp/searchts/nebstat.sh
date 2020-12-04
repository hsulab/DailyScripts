#!/bin/bash
printf "Enter number of images: "
read nimg
for i in `seq 1 $nimg`
do
    # dirname
    dirname=$i
    if [ $i -lt 10 ]
    then
        dirname=0$i
    fi
    # check 
    printf "%s -> " $dirname
    if [ -d $dirname ]
    then
        grep "E0=" $dirname/OSZICAR | tail -n 1
    else
        echo "directroy not exists"
    fi
done
