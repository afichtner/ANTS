#!/bin/bash

# open a file where the names of 0-byte files will be entered. The user can go back and check what happened to that data.

for file in $1/* ;

do
	
	filesize=$(du $file | awk '{print $1}');
	
	if [ $filesize == 0 ] ; then
	
		
		rm $file;
		
	fi;
done
