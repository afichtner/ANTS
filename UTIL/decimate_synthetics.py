import pyasdf
import numpy as np

name_dataset = '/Volumes/cowpox/synthetic.h5'
lowpasscorner = 0.025 # in Hz
decimators = [4, 4, 4]


    
    
#ds = pyasdf.ASDFDataSet(name_dataset)
#print ds
#
#print ds.waveforms['SRC_00000000']
#if lowpasscorner*2 >= ds.waveforms['SRC_00000000'].synthetic[0].stats.sampling_rate / np.product(decimators):
#    msg = 'Not safe: Select lower lowpass corner or aliasing will occur'
#    raise ValueError(msg)
#
#for wf in ds.waveforms:
#    print wf
#    wf.synthetic.taper(type='cosine',max_percentage=0.05)
#    wf.synthetic.filter('lowpass',freq=lowpasscorner,corners=5)
#    for d in decimators:
#        wf.synthetic.decimate(factor=d,no_filter=True)
        
     
def process(st,inv):
    
    for tr in st:
        if lowpasscorner*2 >= tr.\
        stats.sampling_rate / np.product(decimators):
            raise ValueError()
    print st
    st.taper(type='cosine',max_percentage=0.05)
    st.filter('lowpass',freq=lowpasscorner,corners=5)
    for d in decimators:
        st.decimate(factor=d,no_filter=True)
    return st
# Make sure to either use a with statement or delete the reference to the
# data set object at the end. Otherwise it might not be able to properly
# close the file which will stall MPI.
with pyasdf.ASDFDataSet("/Volumes/cowpox/synthetic.h5") as ds:
    ds.process(
#        # Pass the processing function here.
        process_function=process,
        # The output filename. Must not yet exist.
        output_filename="test.h5",
        # Maps the tags. The keys are the tags in the input file and traces
        # with that tag will end up in the output file with the corresponding
        # value.
        # Also note that only data with tags present in this map will be
        # processed. Others will be ignored.
        tag_map={"synthetic": "synth_decimated"},cpu_count=1)

#ds = pyasdf.ASDFDataSet("/Volumes/cowpox/synthetic.h5")
#ds.process(process_function=process,output_filename='test.h5',tag_map={"synthetic": "synth_decimated"},cpu_count=1)