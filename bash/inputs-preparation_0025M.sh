#!/bin/bash

##################################################################
#	AUTOMATIZE INPUT FILES GENERATION V.1.0 (05-20-2019)
# This script is intended to modify in an easy and interactive way
# the electric field of the mdp gromacs file and the topology
# of the CNT providing the value for the el. field and the 
# number of carbon atoms as numbered in the topology of the CNT.
# You can also provide the charge to the carbon atoms by submitting 
# the value of the charge you want to give them before listing
# the carbon atoms as is easily understandable from the outputs.
# N.B NEVERTHELESS INTERACTIVE SELECTION, AS IS INTENDED IN THIS 
# SCRIPT, IS USUALLY AVOIDED SINCE IT HAS NOT MUCH FLEXIBILITY 
# AND CONTROL (I.E. PROVIDING STRINGS INSTEAD OF NUMBERS).
# IN THIS CASE I FOUND IT SIMPLER TO PROVIDE THEM ITERATIVELY 
# INSTEAD OF MAKING A HUGE COMMAND LINE PROVIDING N VALUES THAT 
# COULD EASILY BRING TO ERRORS AND MISTAKES.
# IN CONCLUSION THIS BASH FILE HAS TO BE USED BY PEOPLE WELL 
# AWARE OF THE PROBLEM THIS SCRIPT IS DEALING WITH.
#
#						L. SAGRESTI
#
#
	
set -e

read -p "Please submit the electric field you want provide to the whole system [V/nm] followed by ENTER: " Electric_Field

#Electric_Field=$1

folder="Scripts_0025M/"

OUTPUTMDP="electric_field_try.mdp"

topology_file="topol_ch_7_7_60.itp"

date=$(date +%m-%d-%H-%M)

echo "You are creating an mdp file with "$Electric_Field" V/nm electric field in the z-direction"

carbon_atoms=()

if [ -e "$date""_0025M_elf_""$Electric_Field" ]; then

#	rm -rf "$date""_elf-""$Electric_Field"
	echo "I cannot remove the folder since it was created before with the same name, so pay attention"
	#insert an exit here
fi

mkdir "$date""_0025M_elf_""$Electric_Field" 

cp "$folder"* "$date""_0025M_elf_""$Electric_Field""/"

rm -rf "$date""_1M_elf_""$Electric_Field""/topol_ch_7_7_60_MODIFIED.itp"


sed -e "s/XXXX/$Electric_Field/g" $folder"electric_field_try-XXXX.mdp" > "$date""_0025M_elf_""$Electric_Field""/"$OUTPUTMDP

echo "You are going to modify the CNT topology adding or removing charges at certain C atoms"

read -p "Please submit the total number of carbons you want to change the charge followed by ENTER: " number

read -p "Please submit the total charge you want to give to the previous carbon atoms followed by ENTER: " charge

for i in `seq 1 $number`;
do


	read -p "Please submit the Carbon-ID you want to change followed by ENTER: " first

	carbon_atoms+=($first)

done

cp $folder$topology_file "$date""_0025M_elf_""$Electric_Field""/""tmp_file_0.itp"

n=0
m=1

for i in "${carbon_atoms[@]}"
do


	variable=`awk '{print $1}' $folder$topology_file | grep -n $i | head -1 | tr : " " | awk '{print $1}'`

#	echo $variable

	cd $folder

	A=$(sed -n "${variable}{p;q}" $topology_file)

#	echo "$A"

	B=$(echo "$A" | awk -v charge=$charge '{print "    "$1"         "$2"      "$3"    "$4"      "$5"     "$6"      "charge"     "$8"   "$9,$10,$11}' )

#	echo "$B"

	cd ../"$date""_0025M_elf_""$Electric_Field""/"

	sed "s/$A/$B/g" "tmp_file_"$n".itp" > "tmp_file_"$m".itp"

#	echo "MIAO"

	rm -rf "tmp_file_"$n".itp"

	n=$(expr $n + 1)
	
#	echo "$n"

	m=$(expr $m + 1)

#	echo "$m"	

	cd ../
done

cd "$date""_0025M_elf_""$Electric_Field""/"

mv "tmp_file_"$n".itp" "topol_ch_7_7_60_MODIFIED.itp"


cluster_folder=$(pwd)

sed -e "s#XXXXXX#$cluster_folder#g" submission_grompp_XXXX.sh > submission_grompp.sh

sed -e "s#XXXXXX#$cluster_folder#g" submission_job_XXXX.sh > submission_job.sh

rm -rf *_XXXX.sh

log=$(printf "The $date you have created an mdp file to run 1 ns with an electric field of $Electric_Field V/nm;\nA modified CNT topology charging certain C atoms with charge of $charge named topol_ch_7_7_60_MODIFIED.itp has been created.\nThese inputs could lead to the generation of roughly 150Mb of data after the run.\nGOOD LUCK PAL!!")

echo "$log" > input_preparation.log

#echo "You have successfully created the new modified topology named topol_ch_7_7_60_"$date"_MODIFIED.itp !! GOOD LUCK!!"

cd ../

exit 0


