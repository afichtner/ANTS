#Check the availability of two stations for crosscorrelation
#how...?

def avail_xcorr():
    #checks the files of the current channellist
    #Open the separate files for processing
filelist = open('av_channellist.txt', 'r')
filelist = filelist.readlines()

for i in range(0, len(filelist)-1):
    for j in range(0, len(filelist)-1):
        if i<j:
            #check times of filelist(i) and filelist(j) against one another
            
