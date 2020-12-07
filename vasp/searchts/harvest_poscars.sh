#!/bin/bash
pattern="^[0-9][0-9]"
cwd=`pwd`
tempdir=$cwd/contcars
rm -r $tempdir
mkdir $tempdir
for posdir in `ls`
do
    #echo $posdir | grep "^[0-9].*"
    if [[ $posdir =~ $pattern ]]
    then 
        echo $posdir
        cd $posdir
        pos2cif.pl CONTCAR
        poscarname=CONTCAR_${posdir}.cif
        mv CONTCAR.cif $tempdir/$poscarname
        cp CONTCAR $tempdir/CONTCAR_${posdir}
        cd -
    fi
done
tar -czf temp_contcars.tgz $tempdir
#sz temp_poscars.tgz
