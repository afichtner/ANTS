#!/bin/bash

# open a file where the names of 0-byte files will be entered. The user can fo back and check what happened to that data.
printf 'The following selections returned zero data:\n\n' >> $2

for file in $1/* ;

do
	
	filesize=$(du $file | awk '{print $1}');
	
	if [ $filesize == 0 ] ; then
		
		printf $file >> $2
	
		
		rm $file;
		
	fi;
done
