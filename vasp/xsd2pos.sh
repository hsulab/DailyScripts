#!/bin/bash
#########################################################################
# File Name: xsd2pos.sh
# Author: wdhu@59.78.92.45
# Created Time: Mon 22 Dec 2014 11:54:27 PM CST
# Revised by jyxu at Thomas
#########################################################################
if [ "$1" == "-h" -o "$1" == "--help" ];then
	cat > temp_x2p_help << EOF

[34m#########################################################################
# File Name: xsd2pos.sh
# Author: wdhu
# Created Time: Mon 22 Dec 2014 11:54:27 PM CST
# Revised by jyxu at Thomas
#########################################################################
[0m
Creat inputs for vasp according to the xsdfiles.
The format of POSCAR is vasp5. The creation of POTCAR needs the Pseudopotential Library;
the INCAR vasp.script created or updated need the old files in \$HOME/example or work dir.
the # in INCAR would not replace. 
you can set some parameters by .setting_xsd2pos in home directory or work directory; 
or set from script ARGVs. 
The priority of parameters is 
(1)script ARGVs; (2).setting_xsd2pos in work directory; (3).setting_xsd2pos in home directory.

Examples                               Choices                                                             
[31m cpot=1 [0m                          [34m #0:don't create POTCAR;not 0:create POTCAR [0m 
[31m cscript=1  [0m                      [34m #0:don't create vasp.script;not 0:create; -1 auto [0m 
[31m ckpoint=-1 [0m                      [34m #0:don't create KPOINTS;not 0:create KPOINTS; -1: auto[0m  
[31m setmag=1 [0m                        [34m #0:don't set MAGMOM in INCAR;not 0:set MAGMOM; -1: auto [0m 
[31m setpU=1  [0m                        [34m #0:don't set DFT plus U; not 1:set DFT+U; -1: auto [0m 
[31m sortZ=1  [0m                        [34m #0:sort default ; not 0: sort by Z axis [0m 
[31m posfm=5  [0m                        [34m #the Format of POSCAR: 4 or 5 [0m  
[31m cmark=_c  [0m                        [34m #Constrainmark: _c(default) [0m  
[31m potmode=new [0m                     [34m #new, ori, sv, pv ,GW or Custom  [0m
[31m potdir=\$HOME/pot/potpaw_PBE2010 [0m [34m #Pseudopotential Library directory [0m  
[31m examdir=\$HOME/example [0m           [34m #examples directory [0m 
[31m pUdir=\$HOME/script/Uvalue [0m       [34m #the dir that save U parameters [0m 

The Argv setting
You can set the parameters like the examples above.

[33;1mIf you have any problem or find any BUGs, please contact me.[0m

EOF
	cat temp_x2p_help
	rm -rf temp_x2p_help
	exit 0
fi
rm -rf temp_x2p_*

# Default values for this script 
cpot=1
cscript=1
ckpoint=-1
setmag=1
setpU=1
posfm=5
potmode=new
potdir=$HOME/apps/vasp/pot/potpaw_PBE
examdir=$HOME/example/
pUdir=$HOME/script/Uvalue/
cmark=_c
sortZ=1        #>0:sort by atoms Z axis; 0: only sort by atom element name.

#if [ ! -f *.xsd ];then 
#	rz  #2>&1 > /dev/null
#fi
if [ ! -f *.xsd ];then 
    exit 0
fi


xsdfile=`ls *.xsd | head -1`

echo -e "\n  >>>>>>>>>>>>>>>>>> Check Information <<<<<<<<<<<<<<<<<<\n"

# import the setting files
if [[ -f "$HOME/.setting_xsd2pos" ]];then
	sed -i 's/ //g'  $HOME/.setting_xsd2pos
	. $HOME/.setting_xsd2pos
	echo -e "     Read setting in home directory."
fi
if [[ -f .setting_xsd2pos ]];then
	sed -i 's/ //g' ./.setting_xsd2pos
	. ./.setting_xsd2pos
	echo -e "     Read setting in work directory."
fi

# read ARGV 
until [[ $# -eq 0 ]]
do
	if [[ $1 =~ "=" ]];then
		argvs1=`echo $1 | awk -F = '{print $1}'`
		argvs2=`echo $1 | awk -F = '{print $2}'`
		if [ -z "$argvs1" -o -z "$argvs2" ];then
			echo -e "     Error in Read ARGV, Exit!!!  "
			exit 
		fi
		eval $argvs1=$argvs2	
	elif [[ $1 == "new" ]]; then
		potmode=new
	elif [[ $1 == "ori" ]]; then
		potmode=ori
	elif [[ $1 == "pv" ]]; then
		potmode=pv
	elif [[ $1 == "sv" ]]; then
		potmode=sv
	elif [[ $1 == "GW" ]]; then
		potmode=GW
	elif [[ $1 == "Custom" ]]; then
		potmode=Custom
	elif [[ $1 == "custom" ]]; then
		potmode=Custom
	elif [[ $1 == "cus" ]]; then
		potmode=Custom
	elif [[ $1 == "sortZ" ]]; then
		sortZ=1
	elif [[ $1 == "sortN" ]]; then
		sortZ=0
	fi
	shift
done

if [[  -f "$xsdfile" ]];then
	echo -e "    The xsdfile is : \033[31;1;4m$xsdfile\033[0m \n"
else
	echo -e "    The xsdfile($xsdfile) does not exist!!!"
fi

# create poscar 
elementlist="	H He\
				Li Be  B  C  N  O  F Ne\
				Na Mg Al Si  P  S Cl Ar\
				 K Ca Sc Ti  V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr\
				Rb Sr  Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te  I Xe\
				Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Hf Ta  W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn\
				Fr Ra Ac Th Pa  U Np Pu An Cm Bk Cf Es Fm Md No Rf Db Sg Bh Hs Mt Ds Rg Cn Uut Uuq Uup Uuh Uus Uuo End"

if [[ "`grep Atom3d $xsdfile|grep Components| wc -l`" == 0  ]];then
	echo -e "    None elements information in $xsdfile file. Exit!\n    You can try \"xsd2pos -h\" for help!!!"
fi

grep Atom3d $xsdfile | grep Components | awk -F Name=\" '{print$2}' | awk -F \" '{print $1}' > temp_x2p_Name
grep Atom3d $xsdfile | grep Components | awk -F Components=\" '{print$2}' | awk -F \" '{print $1}'|sed 's/[0-9]//g;s/,//g;s/\.//g;s/+//g' > temp_x2p_Element
grep Atom3d $xsdfile | grep Components | awk -F FormalSpin=\" '{print$2}' | awk -F \" '{print $1+0}'|sed 's/00//g' > temp_x2p_Spin
#grep Atom3d $xsdfile | grep Components | awk -F XYZ=\" '{print$2}' | awk -F \" '{print $1}' |awk -F \, '{printf"%18.12f%18.12f%18.12f\n",$1,$2,$3}' > temp_x2p_XYZ
grep Atom3d $xsdfile | grep Components | awk -F XYZ=\" '{print$2}' | awk -F \" '{print $1}' |awk -F \, '{printf"%18.12f%18.12f%18.12f\n",$1-($1<0?int($1-1):int($1)),$2-($2<0?int($2-1):int($2)),$3-($3<0?int($3-1):int($3))}' > temp_x2p_XYZ
#grep Atom3d $xsdfile | grep Components | awk -F RestrictedProperties=\" '{print$2}' | awk -F \" '{print $1}' | awk '{if($0 ~ "X"){x="F" }else{x="T"};if($0 ~ "Y"){y="F"}else{y="T"};if($0 ~ "Z"){z="F"}else{z="T"};{printf "%4s%4s%4s\n",x,y,z}}' > temp_x2p_TF
grep Atom3d $xsdfile | grep Components | awk -F RestrictedProperties=\" '{print$2}' | awk -F \" '{if($1=="FractionalXYZ"){print "F   F   F"}else{print "T   T   T"}}' > temp_x2p_TF
grep SpaceGroup $xsdfile | awk -F Vector '{printf"%s\n%s\n%s\n",$2,$3,$4}' | awk -F"\"|," '{printf "% 14.8f %14.8f %14.8f\n",$2,$3,$4}' > temp_x2p_lattice

if [[ "$sortZ" -gt 0 ]];then
	paste -d:   temp_x2p_Element  temp_x2p_Name  temp_x2p_XYZ temp_x2p_TF temp_x2p_Spin | awk -F: '{printf"%s:%s\n",  match("'"$elementlist"'",$1),$0}' | sort -n -t: -k1 -k4.37,4.45 -k4.21,4.28| grep -n :  > temp_x2p_tot
else
	paste -d:   temp_x2p_Element  temp_x2p_Name  temp_x2p_XYZ temp_x2p_TF temp_x2p_Spin | awk -F: '{printf"%s:%s\n",  match("'"$elementlist"'",$1),$0}' | sort -n -t: -k1 -s | grep -n :  > temp_x2p_tot
fi


####shortinfo####
awk -F ':' '{print $3}' temp_x2p_tot | uniq -c > temp_x2p_Elenum
awk -F ':' '{printf" %s %s\n", $3,$6}' temp_x2p_tot | awk '{printf" %s  %s%s%s\n",$1,$2,$3,$4}' > temp_x2p_EleTF
rm -f temp_x2p_shortinfo
for i in `cat temp_x2p_Elenum | awk '{print $2}'`
do
	echo  "`grep -w $i temp_x2p_Elenum | awk '{printf "%3s %d\n", $2,$1}'` `grep -w $i temp_x2p_EleTF | grep TTT |wc -l` `grep -w $i temp_x2p_EleTF | grep FFF |wc -l`" >> temp_x2p_shortinfo
done
#######
echo -e "    \033[4m Element       Numbers      \033[0m"
awk 'BEGIN{t=0;f=0;tot=0};{printf"         %2s \033[1;33m%4d\033[0m \033[32m%4d\033[0m(T)\033[31m%4d\033[0m(F)\n",$1,$2,$3,$4}{t+=$3;f+=$4;tot+=$2};END{printf"     ___________________________\n      Total \033[1;33m%4d \033[0m\033[32m%4d\033[0m(T)\033[31m%4d\033[0m(F)\n",tot,t,f}' temp_x2p_shortinfo

####creat POSCAR####
echo -e "AutoCreatByScript\n 1.000000000" > POSCAR
cat temp_x2p_lattice >> POSCAR
awk '{printf "%5s", $1}' temp_x2p_shortinfo >> POSCAR
echo "" >> POSCAR
awk '{print $2}' temp_x2p_shortinfo | awk '{printf "%5d",$1}' >> POSCAR
echo -e "\nSelective Dynamics\nDirect" >> POSCAR
awk -F : '{printf"%s    %s\n", $5, $6}' temp_x2p_tot >> POSCAR
if [[ "$posfm" == 4 ]]; then
	sed 6d POSCAR
fi
#####
###fort.188 Creat###
grep -i "$cmark" temp_x2p_tot | sed 's/:/ /g;s/\t/ /g;s/  / /g'|awk '{printf"%3s %-2s %-5s%10.6f%10.6f%10.6f  %s %s %s\n",$1,$3,$4,$5,$6,$7,$8,$9,$10}' > temp_x2p_Constrain
if [[ -f temp_x2p_Constrain ]];then
	nc=`cat temp_x2p_Constrain | wc -l`
else
	nc=0
fi
if [[ "$nc" -eq 2 ]];then
	awk '{printf "%f %f %f ",$4,$5,$6}' temp_x2p_Constrain  > temp_x2p_ConstrainXYZ
	dis=`$HOME/script/converter/distance temp_x2p_lattice temp_x2p_ConstrainXYZ |awk '{printf"%.5f", $1}'`
	if [[ ! -f "fort.188" ]];then
		echo -e "1\n3\n6\n4\n0.04\n`awk '{printf "%d   ",$1}' temp_x2p_Constrain`  $dis\n0" > "fort.188"
	else
		sed -n '1,5'p fort.188 > temp_x2p_fort188
		echo -e "`awk '{printf "%d   ",$1}' temp_x2p_Constrain`  $dis\n0" >> temp_x2p_fort188
		mv temp_x2p_fort188 fort.188
	fi
	echo -e "\n    find two contrained atoms, shift to TS mode :"
	cat temp_x2p_Constrain | awk '{printf "    %s\n",$0}'
	echo -e "    Its distance : \033[33;1m`echo $dis `\033[0m"
	echo -e "    Creat fort.188 successfully ! \n"
elif [[  "$nc" == 0 ]];then
	echo ""
else
	echo -e "\n    find $nc contrained atom(s), fort.188 don't Creat !\n"
fi
if [[ -f INCAR ]]; then
	if [[ "$nc" -eq 2 ]];then
		sed -i 's/^IBRION *=.*\#/IBRION = 1      #/g' INCAR
		sed -i "s/^NELMIN.*\#/NELMIN = 6        #/g" INCAR
	else
		sed -i 's/^IBRION *=.*\#/IBRION = 2      #/g' INCAR
		sed -i "s/^NELMIN.*\#/NELMIN = 3        #/g" INCAR
	fi
fi

######

####Spin Setting####
if [[ "$setmag" == -1 ]]; then
	if [[ -f INCAR ]]; then
		setmag=1
	else
		setmag=0
	fi
fi
if [[ "$setmag" -gt 0 ]]; then
	if [[ ! -f INCAR ]]; then
		cp $examdir/INCAR .
	fi
	magmomtag=`awk -F : '{print $7 }' temp_x2p_tot | uniq -c | awk '{if($1 == 1){printf" %d",$2}else{if ($2 >= 0 ){printf " %d*%d",$1,$2}else{for(i=1; i<=$1;i++){printf" %d",$2}}}}'`
	if [[ `grep ISPIN INCAR|awk -F \# '{print $1}'|awk -F = '{print$2}'|sed 's/ //g'` == 2 ]];then
		sed -i 's/MAGMOM *=.*/MAGMOM ='"$magmomtag"' #magnetic/' INCAR
		if [[ "`grep MAGMOM INCAR|awk -F \# '{print $1}'`" ]]; then
			echo -e "    The MAGMOM-tag has set !"
			echo -e "    `grep MAGMOM INCAR|awk -F \# '{print $1}'`"
		else
			echo -e "    The MAGMOM-tag do take effect !!"
		fi
	else
		echo -e "    the ISPIN-tag isn't set to 2 !!"
	fi
else
	echo -e "    Don't set MAGMOM-tag !!"
	echo -e "    Spin: $magmomtag .\n"
fi
###################

####Creat POTCAR####
if [[ "$cpot" -gt 0 ]]; then
	if [[ -d "$potdir" ]];then
		rm -f POTCAR temp_x2p_potlog
		if [[ "$potmode" == "new" ]]; then
			for i in `awk '{print $1}' temp_x2p_shortinfo `
			do
				if [[ -f $potdir/"$i"_new/POTCAR_$i ]];then
					cat $potdir/"$i"_new/POTCAR_$i >> POTCAR
					echo ''$i'_new ' >> temp_x2p_potlog
				else
					cat $potdir/"$i"/POTCAR_$i >> POTCAR
					echo "$i " >> temp_x2p_potlog
				fi
			done
		elif [[ "$potmode" == "sv" ]]; then
			for i in `awk '{print $1}' temp_x2p_shortinfo `
			do
				if [[ -f $potdir/"$i"_sv/POTCAR_$i ]];then
					cat $potdir/"$i"_sv/POTCAR_$i >> POTCAR
					echo ''$i'_sv ' >> temp_x2p_potlog
				else
					cat $potdir/"$i"/POTCAR_$i >> POTCAR
					echo "$i " >> temp_x2p_potlog
				fi
			done
		elif [[ "$potmode" == "pv" ]]; then
			for i in `awk '{print $1}' temp_x2p_shortinfo `
			do
				if [[ -f $potdir/"$i"_pv/POTCAR_$i ]];then
					cat $potdir/"$i"_pv/POTCAR_$i >> POTCAR
					echo ''$i'_pv ' >> temp_x2p_potlog
				else
					cat $potdir/"$i"/POTCAR_$i >> POTCAR
					echo "$i " >> temp_x2p_potlog
				fi
			done
		elif [[ "$potmode" == "GW" ]]; then
			for i in `awk '{print $1}' temp_x2p_shortinfo `
			do
				if [[ -f $potdir/"$i"_GW/POTCAR_$i ]];then
					cat $potdir/"$i"_GW/POTCAR_$i >> POTCAR
					echo ''$i'_GW ' >> temp_x2p_potlog
				else
					cat $potdir/"$i"/POTCAR_$i >> POTCAR
					echo "$i " >> temp_x2p_potlog
				fi
			done
		elif [[ "$potmode" == "Custom" ]]; then
			for i in `awk '{print $1}' temp_x2p_shortinfo `
			do
				potname=`grep -w $i $HOME/.potcfg|awk '{print $2}'`
				if [[ -f $potdir/"$potname"/POTCAR_$i ]];then
					cat $potdir/"$potname"/POTCAR_$i >> POTCAR
					echo "$potname " >> temp_x2p_potlog
				else
					cat $potdir/"$i"/POTCAR_$i >> POTCAR
					echo "$i " >> temp_x2p_potlog
				fi
			done
		elif [[ "$potmode" == "ori" ]]; then
			for i in `awk '{print $1}' temp_x2p_shortinfo `
			do
				cat $potdir/"$i"/POTCAR_$i >> POTCAR
				echo "$i " >> temp_x2p_potlog
			done
		fi
		echo -e "    The POTCAR Created! The Elements is : "
		echo -e "    (\033[1;35m`cat temp_x2p_potlog|xargs`\033[0m)\n"
	else
		echo -e "    The POTCAR folder ($potdir) does not exist!!!"
	fi
fi

################
####Creat KPOINTS###
if [[ "$ckpoint" == -1 ]]; then
	if [[ -f KPOINTS ]]; then
		ckpoint=0
	else
		ckpoint=1
	fi
fi
if [[ "$ckpoint" -gt 0 ]]; then
	kpoint=`awk '{print sqrt($1*$1+$2*$2+$3*$3)}' temp_x2p_lattice | awk '{printf"%d ",int(20/$1)+1}'`
	echo -e "mesh auto\n0\nG\n$kpoint\n0 0 0" > KPOINTS
	echo -e "    Creat KPOINTS successfully ( $kpoint) !!!\n"
fi

# Create vasp.script 
if [[ "$cscript" == -1 ]]; then
	if [[ -f vasp.script ]]; then
		cscript=1
	else
		cscript=0
	fi
fi
if [[ "$cscript" -gt 0 ]];then
	if [[ ! -f vasp.script ]]; then
		cp $examdir/vasp.script .
	fi
	sed -i 's/^#PBS -N .*/#PBS -N '"${xsdfile%.*}"'/g' vasp.script
fi

###set DFT+U#######
#$HOME/script/autolda.sh
if [[ "$setpU" == -1 ]]; then
	if [[ -f INCAR ]]; then
		setpU=1
	else
		setpU=0
	fi
fi
if [[ "$setpU" -gt 0 ]]; then
	if [[ ! -f INCAR ]]; then
		cp $examdir/INCAR .
	fi
	if [[ -d "$pUdir" ]]; then
		for i in `awk '{printf"%s \n" , $1}' temp_x2p_shortinfo`
		do
		
			if [[ ! -f "$pUdir/$i" ]]; then
				#echo -e "    None $i U parameters exist, exit!!"
				#echo -e "  <<<<<<<<<<<<<<<<<<<<<<<<=>>>>>>>>>>>>>>>>>>>>>>>>\n"
				#rm -f temp_x2p_* $xsdfile
				#exit 0
				plusUerror=1
				missU=$missU" "$i
			fi
		done
		if [[ "$plusUerror" == 1 ]];then
			echo -e "    The DFT+U parameters of ($missU ) don't exist,\n    DFT+U set cancel! You can try \"xsd2pos -h\" for help!"
		else
			for i in `awk '{printf"%s \n",  $1}' temp_x2p_shortinfo`
			do
				cat $pUdir/$i | awk -F \# '{print $1}'| xargs >>  temp_x2p_Uvalue
			done
			sed -i 's/LDAUL *=.*/LDAUL ='"`cat temp_x2p_Uvalue|awk '{printf" %2s  ",$1}'`"'/'   INCAR
			sed -i 's/LDAUU *=.*/LDAUU ='"`cat temp_x2p_Uvalue|awk '{printf" %3s ",$2}'`"'/'   INCAR
			sed -i 's/LDAUJ *=.*/LDAUJ ='"`cat temp_x2p_Uvalue|awk '{printf" %3s ",$3}'`"'/'   INCAR
			echo -e "    The U parameters updated!!!\n"
		fi
	fi
fi

#END
echo -e "  <<<<<<<<<<<<<<<<<<<<<<<<<<<=>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
rm -f temp_x2p_* $xsdfile 
chmod +x INCAR KPOINTS POSCAR POTCAR vasp.script
