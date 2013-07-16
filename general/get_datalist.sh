#!/bin/bash
#process data in the specified directory
datapath=$1

if [ -a av_daylist.txt ]; then rm av_daylist.txt; fi
if [ -a channellist.txt ]; then rm channellist.txt; fi


ls $datapath | grep 'continuous' > av_daylist.txt
    
    numdays=`cat av_daylist.txt | wc -l`
    
        for i in $(seq 1 $numdays); 
            do
            
            dataday=` awk "NR==$i" av_daylist.txt `
            
            find $datapath/$dataday/BH_RAW/ -type f >> channellist.txt
            find $datapath/$dataday/Resp/ -type f >> resplist.txt
            
            
        done



