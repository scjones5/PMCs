#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.linalg import svd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import socket
import os
import os.path
import sys
import zipfile

print(socket.gethostname())
if socket.gethostname()=="dirac2":
        prefix = "/mnt/tank/home/sjones/PMCs/T-matrix_code/"
        datafile = "/mnt/tank/home/sjones/PMCs/detected_ascii/ss48009/ss48009_90.89.txt"
elif socket.gethostname()=="tank":
        prefix = "/home/sjones/PMCs/T-matrix_code/"
        datafile = "/home/sjones/PMCs/detected_ascii/ss48009/ss48009_90.89.txt"
else: #we must be on rocksteady
        prefix = "/home/sjones/Documents/PMCs/T-matrix_code/"
        datafile = "/home/sjones/Documents/PMCs/detected_ascii/ss48009/ss48009_90.89.txt"

model = "clapp"
#Initial guesses
temp = 135
n_t = 500
shape = 1.9

#here is where we will write the outputs
fileOut = prefix+"occOutputs_011117.txt"
if os.path.isfile(fileOut):
    os.remove(fileOut)

def getPred(parS, observed, hiRes):
    temp = parS[0]
    shape = parS[1]
    n_t = parS[2]
    print("Temperature: "+ str(temp))
    print("EPS parameter: "+ str(shape))
    print("Particle density: "+ str(n_t))
    temp_par = str(round(temp, 1))+"K"
    eps_par = "{0:g}".format(round(shape, 1))+" "
    print("eps_par:", eps_par)
    with open(prefix+"optical_consts/"+model+"/coefficients/calc_coefs_"+temp_par+".txt","r") as n_file:
        wlength = []
        coefArr = []
        eps_line = 6e9 # some arbitrarily large number
        for num, line in enumerate(n_file, 1):
            #print(line)
            if eps_par in line:
                eps_line = num
                continue
            elif num>eps_line:
                if len(line.strip()) == 0:
                    break
                col = line.split()
                #print(col)
                wlength.append(float(col[0]))
                coefArr.append(float(col[1]))
        coefArr2 = np.exp(-(np.asarray(coefArr)*(n_t*(1e-12))*(100000*(1e6))))
    #print("after for loop")
    xx = wlength
    ext = coefArr2
    #print(xx)
    ext_resid = np.interp(hiRes[::-1], xx[::-1], ext)
    xx = np.asarray(xx)
    datafile = prefix+"optical_consts/"+model+"/coefficients/bestFit_spec.txt"
    datafile_id = open(datafile, 'w+')
    bestSpec = np.array([xx, ext])
    bestSpec = bestSpec.T
    np.savetxt(datafile, bestSpec, fmt='%.5f')
    datafile_id.close()
    model_data = (ext_resid)# - offset
    #resid = (abs(savgol_filter(observed, 51, 2) - model_data))
    resid = abs(observed - model_data)/0.1
    wn_loop = 1/(hiRes/10000)
    #plt.clf()
    #plt.plot(wn_loop, savgol_filter(observed, 51, 2))
    #plt.plot(wn_loop, model_data)
    #plt.show()
    #plt.pause(0.0001)
    print(sum(resid)/len(resid))
    return(resid)

ctr = 0
#for filename in os.listdir("/mnt/tank/home/sjones/PMCs/detected_ascii/zips/summer2017"):
with open("/mnt/tank/home/sjones/PMCs/detected_311017.txt", "r") as alts:
    for alt in alts:
        col = alt.split()
        occ = str(col[0])
        height = float(col[1])
        if os.path.isfile("/mnt/tank/home/sjones/PMCs/detected_ascii/summer2017/"+occ[5:]+".zip"):
            with zipfile.ZipFile("/mnt/tank/home/sjones/PMCs/detected_ascii/zips/summer2017/"+occ[5:]+".zip") as myzip:
                names = myzip.namelist()
                print(names)
                for zipfilename in names:
                    #now check if this matches altitude from file
                    hght_toget = round(height, 2)
                    fname = occ[5:]+"_"+hght_toget+".txt"
                    print(fname)
                    with myzip.open(names[names == fname]) as noisy_data:
                        wn = []
                        trans = []
                        for line in noisy_data:
                            col = line.split()
                            wn.append(float(col[0]))
                            trans.append(float(col[1]))

    #filename = 'ss69557.zip'
        filename_par = occ+'.zip' #filename.split('.zip')[0]
        #fname_diff = []
        #for zipfilename in names:
            # find individual file closest to 80km
        #    piece = float((zipfilename.split('_')[1]).split('.txt')[0])
        #    fname_diff.append(np.absolute(piece - 80.0))
        #abs_min = min(np.absolute(fname_diff))
        #print(names)
        #file_min = names[fname_diff.index(abs_min)]

        wn = np.asarray(wn)
        obs_resize = ((wn>=2800) & (wn<=3550))
        wn = wn[obs_resize]
        wlen = np.round((1/(wn*100))*1000000, 5)
        trans = np.asarray(trans)
        trans = trans[obs_resize]
        #if any(x>1.05 for x in trans):
        #    continue
        parStart = [temp,shape,n_t]
        res = least_squares(getPred, parStart, args=(trans, wlen), diff_step=[0.02,0.1,0.1], \
                        bounds=([120,1.1,100],[150,3,800]), verbose=1, xtol=1e-6)
        jac = res.jac
        print("jac:", jac)
        #jacT = np.matrix.transpose(jac)
        #print("jacT:", jacT)
        #innProd = np.dot(jacT, jac)
        #resid_var = np.var(res.fun)
        #print("innProd:", innProd)
        print("variance:", np.var(res.fun))
        error = []
        _, s, VT = svd(res.jac, full_matrices=False)
        threshold = np.finfo(float).eps * max(res.jac.shape) * s[0]
        print("threshold:", threshold)
        print("s:", s)
        print("VT:", VT)
        s = s[s > threshold]
        VT = VT[:s.size]
        pcov = np.dot(VT.T / s**2, VT)
        print("covariance:", pcov)
        error = np.sqrt(np.diag(pcov))
        #if np.linalg.cond(innProd) < 1/sys.float_info.epsilon:
        #    covar = np.linalg.inv(innProd)
        #    covar = resid_var*covar
        #for i in range(len(res.x)):
        #    error.append(np.absolute(pcov[i][i]**0.5))
        #else:
        #    error = [0,0,0]
        print(res.x)
        lineWr = np.asarray((filename_par, res.x[0], error[0], res.x[1], error[1], res.x[2], error[2]), dtype='O')
        print(len(lineWr))
        with open(fileOut, 'ab') as f_handle:
            np.savetxt(f_handle, [lineWr], fmt=("%s", "%.3f", "%.3f", "%.3f", "%.3f", "%.3f", "%.3f"), delimiter=' ', newline=os.linesep)
            f_handle.close()
        datafile = prefix+"optical_consts/"+model+"/coefficients/bestFit_spec.txt"
        bestFit = np.loadtxt(datafile)
        bestCoefs = bestFit[:,1]
        bestWlen = bestFit[:,0]
        wnum = 1/(bestWlen/10000)
        int_bestCoefs = np.interp(wn, wnum, bestCoefs)
        print(wn)
        fig = plt.figure(figsize=(5,4))
        ax = fig.add_subplot(111)
        ax.plot(wn, savgol_filter(trans, 51, 2), 'o')
        ax.plot(wn, int_bestCoefs, '-', linewidth=2)
        fig.savefig("/mnt/tank/home/sjones/PMCs/detected_ascii/Figures/"+filename_par+"_plot.png")
        plt.close(fig)
        #plt.clf()
        #plt.plot(wn, savgol_filter(trans, 51, 2), 'o')
        #plt.plot(wn, (int_bestCoefs), '-', linewidth=2)
        #plt.show()
        ctr = ctr + 1
        #break
