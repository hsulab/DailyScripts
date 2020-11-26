#!/bin/bash
merged_data="MD_merged_data"

echo "# Time  Temp  Ener" > $merged_data

laststep=0
for md_data in $*
do
    echo $md_data "LAST STEP: " $laststep
    awk 'NR!=1 {print $1+'$laststep'" "$2" "$3}' $md_data >> $merged_data
    laststep=`cat $merged_data | tail -n 1 | awk '{print $1}'`
done

# 
laststep=`cat $merged_data | tail -n 1 | awk '{print $1}'`
echo $merged_data "TOT  STEP: " $laststep
echo "Total Steps $laststep in $# MD trajectories."
