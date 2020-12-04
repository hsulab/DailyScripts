#!/bin/bash
if [ x$1 = x ]
then 
    echo "wrong"
    exit 0
fi

eterms=$(grep -A 1 "Step Temp" $1 | tail -n 1)
tnames=('Step' 'Temp' 'E_pair' 'TotEng' 'Press' 'v_eb' 'v_ea'\
    'v_elp' 'v_emol' 'v_ev' 'v_epen' 'v_ecoa' 'v_enb' 'v_et'\
    'v_eco' 'v_ew' 'v_ep' 'v_efi' 'v_eqeq')
#tnames=('Step' 'Temperature' 'E_pair' 'Total Energy' 'Pressure'\
#    'eb bond energy' 'ea atom energy' 'elp lone-pair energy'\
#    'emol moleculer energy(always 0.0)' 'ev valence angle energy'\
#    'epen double-bond valence angle energy' 'ecoa valence angle conjuation energy'\
#    'ehb hydrogen bond energy' 'et torsion energy' 'eco conjugation energy'\
#    'ew Van der Waals energy' 'ep Coulomb energy' 'efi electric field energy'\
#    'eqeq charge equilibrium energy')

count=0
for eterm in $eterms
do
    printf "%30s  %f\n" ${tnames[$count]} $eterm
    ((count ++))
done
