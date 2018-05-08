import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os


######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 14,  
       }  

from_path = './Data_new/'
to_path = './fig_new/' 

######### Parameter you should set ###########
start   =  50  # start time
stop    =  450  # end time
step    =  10  # the interval or step

#  if (os.path.isdir('jpg') == False):
#    os.mkdir('jpg')
######### Script code drawing figure ################
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']

    px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
    py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
    grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
    grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
    work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
    work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
    field_ex = data['Particles/field_ex/subset_high_e/electron'].data/exunit
    field_ey = data['Particles/field_ey/subset_high_e/electron'].data/exunit
    field_bz = data['Particles/field_bz/subset_high_e/electron'].data/bxunit

    px = px [abs(grid_y) < 3.2]
    py = py [abs(grid_y) < 3.2]
    work_x = work_x [abs(grid_y) < 3.2]
    work_y = work_y [abs(grid_y) < 3.2]

    gg = (px**2+py**2+1)**0.5
    theta = np.arctan2(py,px)*180.0/np.pi


#    plt.subplot()
    plt.scatter(theta, gg, c='green', s=5, edgecolors='None', alpha=0.2)

 #   plt.legend(loc='upper right')
    plt.xlim(-180,180)
    plt.ylim(0,400)
    plt.xlabel(r'$\theta\ [degree]$',fontdict=font)
    plt.ylabel(r'$\gamma$',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    #plt.show()
    #lt.figure(figsize=(100,100))
    fig = plt.gcf()
    fig.set_size_inches(8, 6.5)
    fig.savefig(to_path+'angle_gamma'+str(n).zfill(4)+'.png',format='png',dpi=80)
    plt.close("all")

    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
