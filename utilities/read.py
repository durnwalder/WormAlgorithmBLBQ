import numpy as np
import collections
import h5py
import os
import pandas as pd

def get_dataset_keys(f):
    keys = []
    f.visit(lambda key : keys.append(key) if type(f[key]) is h5py._hl.dataset.Dataset else None)
    return keys
def readsim(path, filename, keywords, nkeyword, must_contain=0):
    myDict = collections.defaultdict(list)
    arr = os.listdir(path)
    files = [s for s in arr if filename in s]
    files = [s for s in files if not 'clone' in s]
    files.sort()
    for file in files:
        read = h5py.File(path+file, 'r')
        datasets = get_dataset_keys(read)
        if not must_contain or "simulation/results/Potential_Energy/mean/value" in datasets:
            for keyword in keywords:
                matching = [s for s in datasets if keyword in s]
                matching = [s for s in matching if nkeyword not in s]
                if len(matching) == 1:
                    myDict[keyword].append(read.get(matching[0])[()])
                elif len(matching) > 1:
                    v = [s for s in matching if 'value' in s]
                    e = [s for s in matching if 'mean/error' in s]
                    if len(v) != 0:
                        myDict[keyword].append(read.get(v[0])[()])
                    if len(e) != 0:
                        myDict[keyword+"err"].append(read.get(e[0])[()])
    return myDict
def diff(datal, dataw, obsl, obsw):
    rows_list = []
    for rowl, indl in datal.iterrows():
        for roww, indw in dataw.iterrows():
            if abs(indl['theta_pi']-indw['theta']) < 10e-10 and abs(indl['D_Zeeman']-indw['parameters/D']) < 10e-10:
                rows_list.append({'theta': indw['theta'],'parameters/D': indw['parameters/D'],obsw+'_diff': indw[obsw]-indl[obsl]})
    return pd.DataFrame(rows_list)

