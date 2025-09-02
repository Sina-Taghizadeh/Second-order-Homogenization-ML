#!/bin/bash
SECONDS=0 #for time
mkdir DetailedDir #Dir for all detailed parameters
ConvergenceCriteria=2 #percentage

#for p in 0.4 #poisson ratio
for p in $(seq 0.25 0.03 0.40) #poisson ratio
do
	YangModulus=1000.
	PoissonRatio=$p
	PName=${PoissonRatio:2:5}   #Poisson part of file name
	PoissonFolderName="For${PName}Nu"
	mkdir $PoissonFolderName  #for particular poisson
	
#	for i in 4.82 #strut radius and max of LengthDiscretization
	for i in $(seq 0.2 0.42 4.9) #strut radius and max of LengthDiscretization
	do 
		CubeLength=10.0 #overall length
		StrutRadius=$i
		LengthDiscretization=$StrutRadius #local length for mesh for the fist time
		ForName=${StrutRadius:0:5}   #Radius part of file name
		ForName2=${LengthDiscretization:0:5} #Mesh size part of file name
		FileName="1T${ForName}R${ForName2}M${PName}Nu"
		FinalFileName="For${ForName}R${PName}Nu"
		MyFileName="'$FileName'" #for Salome name
		mkdir $FinalFileName  #for all particular r
		
		#editing salome script
		sed -i '27d' PeriodicMeshScript.py
		sed -i "27i a=$CubeLength" PeriodicMeshScript.py
		sed -i '29d' PeriodicMeshScript.py
		sed -i "29i r=$StrutRadius" PeriodicMeshScript.py
		sed -i '31d' PeriodicMeshScript.py
		sed -i "31i Lw=$LengthDiscretization" PeriodicMeshScript.py
		sed -i '33d' PeriodicMeshScript.py
		sed -i "33i name="$MyFileName"" PeriodicMeshScript.py

		salome -t PeriodicMeshScript.py  # run salome script
		
		#Salome to FEniCS
		MedFile="$FileName.med"
		gmsh -3 -format "msh2" -v 0 $MedFile  #Create msh file with Gmsh 
		dolfin-convert "$FileName.msh" "$FileName.xml" #Convert msh file to final xml files

		python ParaDetermination.py $FileName $CubeLength $FinalFileName $StrutRadius $PoissonRatio $YangModulus $ConvergenceCriteria
		conv=$?
		repeat=1
		
		#removing unessential files
		rm "$FileName.med"
		rm "$FileName.msh"
		rm "$FileName.xml"
		rm "${FileName}_facet_region.xml"
		rm "${FileName}_physical_region.xml"
		
		mv "$FileName.tex" "$FinalFileName"/
		
		while [ $conv -eq 0 ]
		do 
			let "repeat+=1" #simple mathematical operation
			LengthDiscretization=$(echo "(1/$repeat) * $StrutRadius" | bc -l) #more complex mathematical operation and for float numbers
			ForName2=${LengthDiscretization:0:5}
			FileName="${repeat}T${ForName}R${ForName2}M${PName}Nu"
			MyFileName="'$FileName'" #for salome name

			#editing salome script
			sed -i '29d' PeriodicMeshScript.py
			sed -i "29i r=$StrutRadius" PeriodicMeshScript.py
			sed -i '31d' PeriodicMeshScript.py
			sed -i "31i Lw=$LengthDiscretization" PeriodicMeshScript.py
			sed -i '33d' PeriodicMeshScript.py
			sed -i "33i name="$MyFileName"" PeriodicMeshScript.py	

			salome -t PeriodicMeshScript.py
			
			#Salome to FEniCS
			MedFile="$FileName.med"
			gmsh -3 -format "msh2" -v 0 $MedFile   #Create msh file with Gmsh 
			dolfin-convert "$FileName.msh" "$FileName.xml" #Convert msh file to final xml files

			python ParaDetermination.py $FileName $CubeLength $FinalFileName $StrutRadius $PoissonRatio $YangModulus $ConvergenceCriteria
			conv=$?
			
			#removing unessential files
			rm "$FileName.med"
			rm "$FileName.msh"
			rm "$FileName.xml"
			rm "${FileName}_facet_region.xml"
			rm "${FileName}_physical_region.xml"
			
			mv "$FileName.tex" "$FinalFileName"/

		done  #for particular nu and r
		mv "$FileName.tex" "$FinalFileName"/
		mv "$FinalFileName" "$PoissonFolderName"/
	        rm ForCheckingConvergence.txt #After each problem
	done #for particular ne 

	duration=$SECONDS
	echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." > Timefile.txt
	mv Timefile.txt "$PoissonFolderName"/
	mv $PoissonFolderName DetailedDir
done #for all nus and rs
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


