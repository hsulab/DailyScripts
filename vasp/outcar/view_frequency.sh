#!/bin/bash

printout="./print-out"
outcar="./OUTCAR"

if [ -f $printout ]
then
    ### create temporary zpe file
    zpe_file='./zpe_tmp'
    if [ -f ${zpe_file} ]
    then
        echo 'zpe_tmp already existed, thus removed!'
        rm ${zpe_file}
    fi
    freq_outs='./freq_tmp'
    if [ -f ${freq_outs} ]
    then
        echo 'freq_tmp already existed, thus removed!'
        rm ${freq_outs}
    fi

	cal_progress=`grep 'Total' ${printout} | tail -n 1`
    # pro0 remove str total, pro1 now, pro2 target
    pro0=${cal_progress##*Total: }
    pro1=${pro0%%/*}
    pro2=${pro0##*/ }

    if [ ${pro1} -eq ${pro2} ]
    then
        # frequency statistics
        grep 'cm-1' ${outcar} > ${freq_outs}

        grep -v 'f/i' ${freq_outs} | awk '{print "f " $10}' >> ${zpe_file}
        grep 'f/i' ${freq_outs} | awk '{print "i " $9}' >> ${zpe_file}

        nf_all=`cat ${zpe_file} | wc -l`
        nf_real=`cat ${zpe_file} | grep 'f' | wc -l`
        nf_img=`cat ${zpe_file} | grep 'i' | wc -l`
        nf_vibrate=$[ ${nf_all} - 0 ] # for nonlinear molecule

        echo "Stat: ${nf_all} all, ${nf_img} img."

        # calculate ZPE
        cat ${zpe_file} | sort -rn -k 2 | head -n ${nf_vibrate} > ${freq_outs}

        nv_real=`grep 'f' ${freq_outs} | wc -l`
        nv_img=`grep 'i' ${freq_outs} | wc -l`
        ZPE=`grep 'f' ${freq_outs} | awk '{sum+=$2};END {print sum/2000}'`

        echo "Vib: ${nv_real} real, ${nv_img} img."

        rm ${zpe_file}
        rm ${freq_outs}
    else
        freq="${pro1}/${pro2} +-"
        ZPE="np.nan"
    fi
    F_info="${ZPE} eV, ${freq}"
    echo $F_info
else
    echo "There is no print-out!"
fi
