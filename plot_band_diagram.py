from matplotlib import pyplot as plt
import pandas as pd
import simudo.example.fourlayer.sweep_extraction as se
import numpy as np
import os

print(os.getcwd())
os.chdir('example_output/Mar31/a/')
print(os.getcwd())

V = 0.648
#I = 0
I = 0 # 1e-24
#I = 1.6257438673747e-20

con = 0.2
win = 0.03
em = 0.1
absrb = 0.2
bar = 0.1
buf = 0.3

#layers = [con, win,  em,  abs, bar, buf]
#layer_names = ['Contact', Window', 'Emiter', 'Absorber', 'Bar', 'Buffer']
layers = [win,  em,  absrb, bar, buf]
layer_names = ['Window', 'Emiter', 'Absorber', 'Bar', 'Buffer']

dat=se.SweepData('.',f'HMA_V={V}')
adict = dict(file=f'HMA_V={V}.plot_meta.yaml')
string = f"V={V}" # will be used in filename of .png and in title of plot. It will be saved to directory the file in SweeData() is in

#dat=se.SweepData('.',f'HMA_I={I}')
#adict = dict(file=f'HMA_I={I}.plot_meta.yaml')
#string = f"I={I}" # will be used in filename of .png and in title of plot. It will be saved to directory the file in SweeData() is in

save_name = f"Bands_"+string+".png"
spatial = dat.get_spatial_data(adict)

total = 0
plt.title(f'HMA band diagram at ' + string)
for l in layers:
    total+=l
    plt.vlines(total, -2, 2.5, linestyle = '--')
plt.plot(spatial["coord_x"], spatial["Ephi_CB"], color="k", label = 'CB')
plt.plot(spatial["coord_x"], spatial["Ephi_VB"], color="k", label = 'VB')
plt.plot(spatial["coord_x"], spatial["Ephi_IB"], color="k", label = 'IB')
plt.plot(spatial["coord_x"], spatial["qfl_CB"], color="b", label = 'CB qfl')
plt.plot(spatial["coord_x"], spatial["qfl_VB"], color="r", label = 'VB qfl')
plt.plot(spatial["coord_x"], spatial["qfl_IB"], color="g", label = 'IB qfl')
plt.legend()
plt.xlabel(r"x ($\mu$m)")
plt.ylabel(r"Energy (eV)")
#plt.ylim([1e135,1e137])
plt.grid()
plt.savefig(save_name)
plt.show()
#"""

##########################################################################
# Plotting the difference in qfl
"""
plt.figure()
plt.title(f'HMA band diagram at ' + string)
total = 0
for l in layers:
    total+=l
    plt.vlines(total, -2, 2.5, linestyle = '--')
dif_qfl = [i-j for i,j in zip(spatial["qfl_CB"], spatial["qfl_VB"]) ]
plt.plot(spatial["coord_x"], dif_qfl, color="orange", label = 'CB_qfl - VB_qfl')
plt.plot(spatial["coord_x"], spatial["qfl_CB"], color="b", label = 'CB qfl')
plt.plot(spatial["coord_x"], spatial["qfl_VB"], color="r", label = 'VB qfl')
plt.plot(spatial["coord_x"], spatial["qfl_IB"], color="g", label = 'IB qfl')
plt.legend()
plt.xlabel(r"x ($\mu$m)")
plt.ylabel(r"Energy (eV)")
#plt.ylim([1e135,1e137])
plt.grid()
#plt.savefig(save_name)
#plt.show()
"""
##########################################################################
# Plotting the difference in band edges
plt.figure()
plt.title(f'HMA band diagram at ' + string)
total = 0
for l in layers:
    total+=l
    plt.vlines(total, -2, 2.5, linestyle = '--')
dif_edges = [i-j for i,j in zip(spatial["Ephi_CB"], spatial["Ephi_VB"]) ]
plt.plot(spatial["coord_x"], dif_edges, color="orange", label = 'E_VB - E_CB')
plt.plot(spatial["coord_x"], spatial["Ephi_CB"], color="k", label = 'E_CB')
plt.plot(spatial["coord_x"], spatial["Ephi_VB"], color="k", label = 'E_VB')
plt.legend()
plt.xlabel(r"x ($\mu$m)")
plt.ylabel(r"Energy (eV)")
plt.grid()


total = 0
dif_edges = np.asarray(dif_edges)
for i in range(len(layers)):
    print('\nAveraged bandgap energy in', layer_names[i])

    a = spatial["coord_x"][spatial["coord_x"] > total]

    e = dif_edges[list(spatial["coord_x"] > total)]
    e = e[a < total+layers[i]]

    a = a[a < total+layers[i]]

    print(round(sum(e)/len(e), 4), 'eV')
    total += layers[i]

plt.show()
#"""