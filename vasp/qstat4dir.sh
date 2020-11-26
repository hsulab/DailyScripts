#!/bin/bash
cur_dir=`pwd`
ifconverg=" reached required accuracy - stopping structural energy minimisation"

for dir in `ls $cur_dir`
do  
    printf "%-20s" "$dir"

    cur_printout=$cur_dir/$dir/print-out
    if [ -f $cur_printout ]
    then 
        lastline=`cat $cur_printout | tail -n 1` 
        if [ "$lastline" = "$ifconverg" ] 
        then 
            echo -ne "Has converged."
        else 
            echo -ne "Not converged."
        fi
    else 
        echo -ne "No print-out yet."
    fi

    calc_number=`ls $dir | grep "\.o" | wc -l`
    echo " Calculated $calc_number times."
done
