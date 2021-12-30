import numpy as np
import matplotlib.pyplot as plt
import h5py

numQuadrants = 4

# trans_spec_data = h5py.File('trans_spec_data.mat','r')
tsd = h5py.File('trans_spec_data.mat','r')
wlVals = np.array(tsd['E_fm0']['lambda'])[0]

def get_monitors(monitor_name,monitor_var,num_monitors):
    # Say you have 4 monitors fm_0,fm_1,fm_2,fm_3. Gets all their data (wrt the relevant variable) and stores it in an
    # array of length 4 called fm
    monitor_data = []

    for i in np.arange(0,num_monitors):
        data = np.array(tsd.get(monitor_name+str(i)).get(monitor_var))
        data = data['real'] + np.complex(0, 1) * data['imag']
        monitor_data.append(data)
    monitor_data = np.asarray(monitor_data)

    return monitor_data

fm = get_monitors('E_fm', 'E', 4)
tm = get_monitors('E_tm', 'E', 4)
fp = get_monitors('E_fp', 'E', 1)