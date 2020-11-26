#!/bin/bash
merged_data="merged_out.ani"

laststep=0
for outani in $*
do
    echo $outani "LAST STEP: " $laststep
    awk -v s="$laststep" '{if (match($1, "STEP")) print $1" = "$3+s; else print $0}' $outani >> $merged_data
    laststep=`grep "STEP" $merged_data | tail -n 1`
    laststep=${laststep#*=}
done

laststep=`grep "STEP" $merged_data | tail -n 1`
echo "Total Steps $laststep in $# MD trajectories."
