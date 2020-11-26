#########################################################################
# File Name: gDns.sh
# Author: jyxu
# mail: ahcigar@foxmail.com
# Created Time: Fri 01 Mar 2019 09:45:45 AM CST
#########################################################################

#!/bin/bash

# check tmp file
cur_dir=`pwd`
qstat_tmp=${cur_dir}/qstat_tmp

if [ -f $qstat_tmp ]
then
    rm $qstat_tmp
    echo "qstat_tmp already existed"
fi

# get qstat info
nstat=`qstat | wc -l`
ntask=$[nstat-2]
qstat | tail -n $ntask > $qstat_tmp

tasks=(`awk '{print $1}' $qstat_tmp`)
for task in ${tasks[@]}
do
    stat=`grep $task $qstat_tmp`
    if [ "`echo $stat | awk '{print $5}'`" = "R" ]
    then
        work_dir=`qstat -f $task | grep 'init' | awk '{print $3}'`
        echo ${work_dir}
        cd $work_dir
        gDn | tail -n 2
    fi
done

# back dir and rm tmp file
cd $cur_dir
rm $qstat_tmp
