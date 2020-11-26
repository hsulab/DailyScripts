#!/bin/sh
###
echo "This Script will delete ALL TASKS after the TaskNumber."
###
if [ $# -eq 0 ]
then
    echo "$0 [TaskNumber]"
    exit
elif [ $# -eq 1 ]
then
    TaskNumber=$1
    qstat | grep 'master' | awk '{print $1}' >> qstat_temp
    for line in `cat qstat_temp`
    do
	    number=${line%.*}
	    if [ $((number)) -gt $((TaskNumber)) ]
	    then
		    qdel $number
	    fi
    done
    rm qstat_temp
else
    echo "$0 [TaskNumber]"
    exit
fi
###
