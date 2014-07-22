from __future__ import print_function
from obspy import Stream


def rotate_streams(str1,str2,baz,rot='NE->RT',verbose=False,outfile=None):
    
    if len(str1) == 0 or len(str2) == 0:
        msg = 'One or both streams are empty. Cannot rotate'
        raise ValueError(msg)
    
    if str1[0].id[-1] in ['N','E'] == False or str2[0].id[-1] in ['N','E'] == False:
        msg = 'Wrong component: Can only rotate \'NE\'.'
        raise ValueError(msg)
        
    
    (str1new,str2new,perck) = find_common_segments(str1,str2)
    
    if perck > 75 and verbose:
        print('Warning! Less than 75 % of original streams could be rotated.',file=outfile)
    
    str = str1new + str2new
    str.rotate(rot,baz)

    return str


def find_common_segments(str1,str2,verbose=False):
    
    if len(str1) == 0 or len(str2) == 0:
        msg = 'One or both streams are empty.'
        raise ValueError(msg)
    
    str1new = Stream()
    str2new = Stream()
    
    n1 = 0
    n2 = 0
    
    numsamp1 = 0
    numsampnew = 0
    
    for i in range(len(str1)):
        numsamp1 += len(str1[i].data)
    
    while n1 < len(str1) and n2 < len(str2):
        
        start1 = str1[n1].stats.starttime
        start2 = str2[n2].stats.starttime
        
        end1 = str1[n1].stats.endtime
        end2 = str2[n2].stats.endtime
        
        # Test if segments overlap at all
        
        if start1 > end2:
            n2 += 1
            continue
        if start2 > end1:
            n1 += 1
            continue
        
        # Find start and endtime
        
        start = max(start1,start2)
        end = min(end1, end2)
        
        str1new += str1.slice(starttime=start,endtime=end)
        str2new += str2.slice(starttime=start,endtime=end)
        
        # increase index
        
        if end1 > end2:
            n2 += 1
        elif end1 == end2:
            n1 += 1
            n2 += 1
        else:
            n1 += 1
        
    for i in range(len(str1new)):
        numsampnew += len(str1new[i].data)
    percentkept = numsamp1/numsampnew * 100 
    
    return(str1new,str2new,percentkept)
    


