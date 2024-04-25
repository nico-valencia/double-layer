# import pandas as pd 
# import numpy as np 
# import matplotlib.pyplot as plt 
# # import scipy.stats as scistats
# # import scipy.optimize as sci 
# # import cartopy.crs as ccrs
import splitwavepy as sw
import numpy as np
# import pandas as pd 
import matplotlib.pyplot as plt
# # import pyvista as pv 
# # import gc 
# # import math 
# import copy 
# from obspy.signal.filter import envelope
# from obspy import read 
# from scipy.interpolate import griddata
# import glob
# from obspy import geodetics
# import itertools 
# from multiprocessing import Pool 
# import multiprocessing
# import os 
from support import plot_trace  #plot_eigen_measure, plot_trans_measure , lambda_filter_2 , plotdata
# import random 



# baz = 230

def baztopolar(baz):
    polar = 0 
    polar = (90-baz) % 360 
    return polar

ang1 = -90 
ang2 = 90 
dt1 = 0 
dt2 = 2

# baz = random.randint(0,360)
baz = 90
print(baz)
# pair    = sw.Pair(split=[(60,0.5),(60.0,0.8)], noise=0.005, pol=baztopolar(baz), delta=0.05)

# wavelet = sw.Pair(split=[(30,0.5),(0.0,0.0)], noise=0.005, pol=baztopolar(baz), delta=0.05)
wavelet = sw.Pair(split=[(40,0.9),(0.0,0.0)], noise=0.01, pol=baz, delta=0.01,width=150,nsamps=8001)
wavelet.rotateto(baz)
# wavelet.unsplit(40,0.9)

kwargs = {'figname': 'synth_trace.png'}
plot_trace(wavelet,baz, **kwargs)
plt.close()
n = len(wavelet.x)
d = 0.01
t = np.linspace(0.0, n*d, n, endpoint=False)
plt.plot(t,wavelet.x,label='x',color='blue')
plt.plot(t,wavelet.y,label='y',color='red')
plt.legend()
plt.show()
plt.close()
# wavelet.unsplit(30,0.5)
# wavelet.unsplit()
# kwargs = {'figname': 'synth_trace.png'}
# plot_trace(wavelet,baz, **kwargs)



# filename = './synthetic_data_fr.txt'

# file = open(filename, "w")
# file.write("name"           +"|"+ 
#            "#Network"       +"|"+ 
#            "Station"        +"|"+ 
#            "Latitude"       +"|"+ 
#            "Longitude"      +"|"+ 
#            "fast1"          +"|"+ 
#            "lag1"           +"|"+ 
#            "fast2"          +"|"+ 
#            "lag2"           +"|"+ 
#            "dist"           +"|"+ 
#            "az"             +"|"+ 
#            "baz"            +"|"+
#            "event"          +"|"+
#            "evlat"          +"|"+
#            "evlon"          +"|"+
#            "evdepth"        +"|"+
#            "evmag"          +"|"+
#            "SNR"            +"|"+
#            "Emin"           +"|"+
#         #    "dfast"          +"|"+
#         #    "dlag"           +"|"+
#            "radsnr"         +"|"+
#            "transsnr"       +"|"+
#            "radenv"         +"|"+   
#            "transenv"       +"|"+ 
#            "\n")

# bazrange = np.arange(0,360,10)
# # for i in range(0,len(path_list),1):
# for i in range(0,len(bazrange),1):
#     print(bazrange[i])
#     # try:
#     if 1 ==1:  

#         wavelet = sw.Pair(split=[(30,0.5),(0.0,0.0)], noise=0.001, pol=bazrange[i], delta=0.05)
#         x = wavelet.splitting_intensity()
#         # print(len(wavelet.x))
#         wave = copy.deepcopy(wavelet)
#         ogwave = copy.deepcopy(wavelet)
#         R = ogwave.x
#         T = ogwave.y
#         Rmod = wave.x
#         Tmod = wave.y
#         print(x)





#         ####### doing the splitting measurement #######

#         #####################

        
#         # define parameter space 

#         param_grid = {
#             'angle1': np.around(np.linspace(ang1,ang2,31),3),
#             'delay1': np.around(np.linspace(dt1,dt2,21),3),
#             'angle2': np.around(np.linspace(ang1,ang2,31),3),
#             'delay2': np.around(np.linspace(dt1,dt2,21),3),
#                     }

#         # print(param_grid)

#         param_combination = list(itertools.product(*param_grid.values()))


#         ###############
#         def valfinder(p):
#             global Rmod,Tmod,R,T
#             ang1 = p[0]
#             dt1 = p[1]
#             ang2 = p[2]
#             dt2 = p[3]


#             # wavelet = pair
#             wave = copy.deepcopy(wavelet)
#             ogwave = copy.deepcopy(wavelet)
            
#             wave.rotateto(bazrange[i])
#             ogwave.rotateto(bazrange[i])  #should this be rotated as well or not? 

#             wave.unsplit(ang1,dt1)
#             wave.unsplit(ang2,dt2)


#             ######## energy maximization method (valencia) ########
#             rad = wave.x
#             trans = wave.y
#             trans_og = ogwave.y
#             rad_energy = np.sum(np.square(rad))
#             trans_energy = np.sum(np.square(trans_og[0:len(trans)] - trans[0:len(trans)]))
#             E = np.sqrt(trans_energy/rad_energy)
#             #######################################################

#             # ######## modified energy minimization method ###########
#             # R = ogwave.x
#             # T = ogwave.y

#             # Rmod = wave.x
#             # Tmod = wave.y 

#             # R_energy = np.sum(np.square(R))
#             # Tmod_energy = np.sum(np.square(Tmod))
#             # E = np.sqrt(Tmod_energy / R_energy )
#             # ##########################################################





#             # E = 0 
#             row = [np.around(ang1,4),np.around(dt1,4),np.around(ang2,4),np.around(dt2,4),np.around((np.max(E)),6)]
#             # print(row)
            
#             return E

#         ######### SNR (Restivo & Helffrich) #########
#         signal = np.max(np.abs(Rmod))
#         noise = 2 * np.std(Tmod)
#         snr = signal/noise
#         #############################################


#         ########## signal/noise constraints ##########
#         # radial signal 
#         rad_noise               = R[0:300]            # 0 - 50 seconds         
#         rad_signal              = R[320:800]         # 54 - 106 seconds 
#         maxradnoise             = np.max(np.abs(rad_noise))
#         maxradsignal            = np.max(np.abs(rad_signal))
#         rad_signal_to_noise     = maxradsignal / maxradnoise 

#         # rad_snr_str             = str(rad_signal_to_noise)

#         # transversal signal 
#         trans_noise             = T[0:300]         # 0 - 50 seconds 
#         trans_signal            = T[320:800]      # 54 - 106 seconds 
#         maxtransnoise           = np.max(np.abs(trans_noise))
#         maxtranssignal          = np.max(np.abs(trans_signal))
#         trans_signal_to_noise   = maxtranssignal / maxtransnoise 

#         # trans_snr_str           = str(trans_signal_to_noise)
        
#         ### envelope filtering ###

#         # radial envelope
#         rad_noise_envelope      = envelope(rad_noise)
#         rad_signal_envelope     = envelope(rad_signal)
#         rad_noise_chop          = rad_noise_envelope[0:980]            # originally has 0 - 1000 samples 
#         rad_signal_chop         = rad_signal_envelope[20:1020]          # originally has 1080 - 2120 samples (1,040 total)
#         rad_env                 = np.max(rad_signal_chop) / np.max(rad_noise_chop)
#         # rad_env_str             = str(rad_env)


#         # transversal envelope
#         trans_noise_envelope    = envelope(trans_noise)
#         trans_signal_envelope   = envelope(trans_signal)
#         trans_noise_chop        = trans_noise_envelope[20:980]
#         trans_signal_chop       = trans_signal_envelope[20:1020]
#         trans_env               = np.max(trans_signal_chop) / np.max(trans_noise_chop)
#         # trans_env_str           = str(trans_env)



#         ##############################################

#         def parallel_grid_search(param_combination):
#             return param_combination, valfinder(param_combination)


#         if __name__ == '__main__':
#             # Use multiprocessing Pool to distribute tasks across multiple processes
#             with Pool(processes=64) as pool:
#                 # Perform parallel grid search
#                 results = pool.map(parallel_grid_search, param_combination)

            
#             # Aggregate results
#                 best_params, best_performance = None, float('0')
#                 data = np.zeros([len(results),5])

#                 for i in range(0,len(results)):
#                     params , performance = results[i] 

#                     if performance > best_performance:
#                         best_params, best_performance = params, performance

#                     data[i] = params[0],params[1],params[2],params[3],performance

                
#             sorted_data = np.argsort(data[:,4])
#             print(data[sorted_data])

#             df = pd.DataFrame(data)
#             df.to_csv('method2.csv',header=None, index=None)
                

#             # Print or return the best parameters and performance
#             print("Best parameters:", best_params)
#             print("Best performance:", best_performance)
#             print(len(results))
#             print(len(data))

#         file = open(filename, "a")
#         file.write(
#                 str(best_params[0]) +"|"+
#                 str(best_params[1]) +"|"+
#                 str(best_params[2]) +"|"+
#                 str(best_params[3]) +"|"+
                

#                 str(snr)  +"|"+

#                 str(best_performance) +"|"+ 
#                 str(rad_signal_to_noise) +"|"+
#                 str(trans_signal_to_noise) +"|"+
#                 str(rad_env) +"|"+
#                 str(trans_env)
#                     + "\n")





#                 #    str(best_params) +"|"+ 
#                 #    str(best_performance) +"|"+ 
#                 #    str(baz) + "\n")
                    
#         print("Saving Measurement to File:",filename)
#     # except:
#     #      print("measurement did not work")

