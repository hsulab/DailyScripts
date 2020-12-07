#!/bin/bash
pattern="^[0-9][0-9]"
cwd=`pwd`
tempdir=$cwd/poscars
rm -r $tempdir
mkdir $tempdir
for posdir in `ls`
do
    #echo $posdir | grep "^[0-9].*"
    if [[ $posdir =~ $pattern ]]
    then 
        echo $posdir
        cd $posdir
        pos2cif.pl POSCAR
        poscarname=POSCAR_${posdir}.cif
        mv POSCAR.cif $tempdir/$poscarname
        cd -
    fi
done
tar -czf temp_poscars.tgz $tempdir
sz temp_poscars.tgz
