#!/bin/bash
#process data in the specified directory
resppath=$1
resploc=$2

if [ -a av_daylist.txt ]; then rm av_daylist.txt; fi
if [ -a resplist.txt ]; then rm resplist.txt; fi


ls $resppath | grep 'continuous' > av_daylist.txt
    
    numdays=`cat av_daylist.txt | wc -l`
    
        for i in $(seq 1 $numdays); 
            do
            
            respday=` awk "NR==$i" av_daylist.txt `
    
            find $resppath/$respday/Resp/ -type f >> resplist.txt
            
            
        done


    numresp=`cat resplist.txt | wc -l`
        
        for i in $(seq 1 $numresp);
        do
        respfile=` awk "NR==$i" resplist.txt `
        starttime=`cat $respfile | grep 'Start date'`
        year=`echo $starttime | grep -o '[0-9]\{4\}'`
        
        respdirname="Resp_$year"
        
        if [ ! -d "$resploc" ]; then
        mkdir $resploc
        fi
        
        if [ ! -d "$resploc/$respdirname" ]; then
        mkdir $resploc/$respdirname
        fi
                
        cp $respfile $resploc/$respdirname
        done
