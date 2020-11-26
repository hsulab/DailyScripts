#!/bin/bash

# check tmp file
cur_dir=`pwd`
echo "Current workding dir is $cur_dir !"

# dirs 
for dir in `ls`
do 
    ts_dir=$dir/ts
    tsra_dir=$dir/ts_ra
    if [ -d $dir ]
    then 
        # ts dir
        echo -n "$ts_dir     " 
        if [ -d $ts_dir ]
        then 
            if [ -f $ts_dir/OUTCAR ]
            then 
                cd $ts_dir
                info=`gDn | tail -n 2`
                echo $info | awk '{print $9 $10 $11}' 
                cd $cur_dir
            else 
                echo "Not calculated."
            fi
        else 
            echo "No dir."
        fi
        # tsra dir 
        echo -n "$tsra_dir     " 
        if [ -d $tsra_dir ]
        then 
            if [ -f $tsra_dir/OUTCAR ]
            then 
                cd $tsra_dir
                info=`gDn | tail -n 2`
                echo $info | awk '{print $9 $10 $11}' 
                cd $cur_dir
            else 
                echo "Not calculated."
            fi
        else 
            echo "No dir."
        fi
    fi
done 
