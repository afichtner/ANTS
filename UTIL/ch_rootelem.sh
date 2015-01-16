#!/bin/bash
sed '
/^<FDSNStationXML/ c\
	<FDSNStationXML>' $1 > temp.txt

mv temp.txt $1

sed '
/^<\/FDSNStationXML/ c\
	<\/FDSNStationXML>' $1 > temp.txt

sed  's/iris\:/iris/g' $1 > temp.txt


mv temp.txt $1


