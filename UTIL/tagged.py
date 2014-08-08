import fnmatch
import os
import sys

from glob import glob


if __name__=='__main__':
    import tagged as t
    todo = str(sys.argv[1])
    tag = str(sys.argv[2])
    
    if todo == '-d':
        t.delete_tagged(tag)
    elif todo == '-mv':
        t.move_tagged(tag,str(sys.argv[3]))
    
    

def delete_tagged(tag):
    """
    Delete all data that carries a certain tag, with all related metadata. 
    Everything carrying that tag will be erased!
    
    """

    # find data
    data = finddata(tag)
    
    
    print '\n\n================== Files to delete: =======================================\n'
    print data
   
    dec = raw_input('\n==================Are you sure [yes/no]?================================\n')
    if dec == 'yes':
        for file in data:
            os.system('rm '+file)
    
    
def move_tagged(tag, destination):
    
    if destination[-1]!='/':
        destination+='/'
    
    data = finddata(tag)
    if os.path.exists(destination+tag) == False:
        os.mkdir(destination+tag)
        
    print '\n\n================== Files to move: =======================================\n'
    print data
   
    dec = raw_input('\n==================Are you sure [yes/no]?================================\n')
    if dec == 'yes':
        for file in data:
            if os.path.exists(destination+tag+'/'+file) == False:
                os.system('mv ' + file + ' ' + destination+tag+'/')
    
def finddata(tag):
    
    data = []
    for root, dirnames, filenames in os.walk('./DATA/'):
        for filename in fnmatch.filter(filenames, '*.'+tag+'.*'):
            data.append(os.path.join(root, filename))
    return data
    
