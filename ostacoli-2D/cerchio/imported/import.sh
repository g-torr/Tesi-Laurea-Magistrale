 #!/bin/bash 
FILE="path"
for nomefile in $(cut -d: -f1 $FILE);
do	echo $nomefile
	nomefile2=$(echo $nomefile|sed -r 's/.inside.cu$//g')
	hostpath=$nomefile2
	echo $hostpath
	mkdir $hostpath
	scp -r $echo'roberto@dev.phys.uniroma1.it:/home/roberto/Documents/Giuseppe/cerchio/changin_R_HigherL/'$nomefile  $hostpath;
done

