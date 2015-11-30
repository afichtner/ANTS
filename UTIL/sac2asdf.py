import pyasdf
from glob import glob
from obspy import read, UTCDateTime

data=pyasdf.ASDFDataSet(filename='OUTPUT_FILES/synthetics_z.h5')

traces = glob('OUTPUT_FILES/SRC.*.MXZ.sem.sac')

for trace in traces:
    tr = read(trace)
    tg = trace.split('/')[-1]
    tg = tg.split('.')[1]
    tr[0].stats.starttime=UTCDateTime(2000,01,01)
    data.add_waveforms(tr,tag=tg)
