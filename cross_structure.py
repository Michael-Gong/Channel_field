#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
  
if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     1.0e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('red')
  c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
  c_orange = matplotlib.colors.colorConverter.to_rgba('orange')
  c_green = matplotlib.colors.colorConverter.to_rgba('green')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128) 
  cmap_go = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_green,c_white_trans,c_orange],128) 
##end for transparent colorbar##
 
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####

  ######### Parameter you should set ###########
  start   =  20  # start time
  stop    =  20  # end time
  step    =  1  # the interval or step
  ######### Script code drawing figure ################

  dict_list = {'00':1999,'01':1875,'02':1825,'05':1725,'10':1600,'20':1400}

  for n in range(start,stop+step,step):
        ywidth=80
         
        plt.subplot(2,3,1)
        #### header data ####
        data = sdf.read('./Data_n/'+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        x  = data['Grid/Grid_mid'].data[0]/1.0e-6
        print('ok')
        y  = data['Grid/Grid_mid'].data[1]/1.0e-6
        X, Y = np.meshgrid(x, y) 
        name = '00'
        den_out = data['Derived/Number_Density/electron'].data/denunit
        den_out = den_out[250:dict_list[name],300-ywidth:300+ywidth]
        den_in = np.zeros_like(den_out)
        ex = (data['Electric Field/Ey'].data-data['Electric Field/Ey_averaged'].data)/exunit
        ex = ex[250:dict_list[name],300-ywidth:300+ywidth]
        ex2 = ex**2

        plt.plot(y[300-ywidth:300+ywidth], den_out.sum(axis=0)/den_out[:,0].size, linewidth=3,linestyle='-',color='black', label=r'$n_{out}$')
        plt.plot(y[300-ywidth:300+ywidth], den_in.sum(axis=0)/den_in[:,0].size, linewidth=3,linestyle='-',color='red', label=r'$n_{in}$')
        plt.plot(y[300-ywidth:300+ywidth], 0.5*(den_in.sum(axis=0)/den_in[:,0].size+den_out.sum(axis=0)/den_out[:,0].size), linewidth=3,linestyle='--',color='orange', label=r'$\frac{n_{out}+n_{in}}{2}$')
        plt.plot(y[300-ywidth:300+ywidth], np.zeros_like(y[300-ywidth:300+ywidth])+0.0/2, linewidth=3,linestyle='--',color='green', label=r'$\frac{n_{ion}}{2}$')
        plt.legend(loc='upper left', fontsize=20, framealpha=0.1)
        plt.xlabel('Y [$\lambda_0$]',fontdict=font)
        plt.ylabel(r'$n_e\ [n_c]$',fontdict=font, color='black')
        plt.ylim(0,1.1)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

        plt1 = plt.twinx()
        plt1.plot(y[300-ywidth:300+ywidth], (ex2.sum(axis=0)/ex2[:,0].size)**0.5 ,linewidth=3, linestyle='-',color='blue')
        plt1.set_ylabel(r'$\sqrt{\bar{E_y^2}}\ [m_ec\omega_0/e]$', color='blue')
        plt1.tick_params('y', colors='blue')
        plt1.set_ylim((0,15))
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);
        

        plt.subplot(2,3,2)
        #### header data ####
        data = sdf.read('./Data_n01/'+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        x  = data['Grid/Grid_mid'].data[0]/1.0e-6
        print('ok')
        y  = data['Grid/Grid_mid'].data[1]/1.0e-6
        X, Y = np.meshgrid(x, y) 
        name = '01'
        den_out = data['Derived/Number_Density/electron'].data/denunit
        den_out = den_out[250:dict_list[name],300-ywidth:300+ywidth]
        den_in = data['Derived/Number_Density/electronin'].data/denunit
        den_in = den_in[250:dict_list[name],300-ywidth:300+ywidth]
        ex = (data['Electric Field/Ey'].data-data['Electric Field/Ey_averaged'].data)/exunit
        ex = ex[250:dict_list[name],300-ywidth:300+ywidth]
        ex2 = ex**2

        plt.plot(y[300-ywidth:300+ywidth], den_out.sum(axis=0)/den_out[:,0].size, linewidth=3,linestyle='-',color='black', label=r'$n_{out}$')
        plt.plot(y[300-ywidth:300+ywidth], den_in.sum(axis=0)/den_in[:,0].size, linewidth=3,linestyle='-',color='red', label=r'$n_{in}$')
        plt.plot(y[300-ywidth:300+ywidth], 0.5*(den_in.sum(axis=0)/den_in[:,0].size+den_out.sum(axis=0)/den_out[:,0].size), linewidth=3,linestyle='--',color='orange', label=r'$\frac{n_{out}+n_{in}}{2}$')
        plt.plot(y[300-ywidth:300+ywidth], np.zeros_like(y[300-ywidth:300+ywidth])+0.1/2, linewidth=3,linestyle='--',color='green', label=r'$\frac{n_{ion}}{2}$')
        plt.xlabel('Y [$\lambda_0$]',fontdict=font)
        plt.ylabel(r'$n_e\ [n_c]$',fontdict=font, color='black')
        plt.ylim(0,1.1)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

        plt1 = plt.twinx()
        plt1.plot(y[300-ywidth:300+ywidth], (ex2.sum(axis=0)/ex2[:,0].size)**0.5 ,linewidth=3, linestyle='-',color='blue')
        plt1.set_ylabel(r'$\sqrt{\bar{E_y^2}}\ [m_ec\omega_0/e]$', color='blue')
        plt1.tick_params('y', colors='blue')
        plt1.set_ylim((0,15))
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);


        plt.subplot(2,3,3)
        #### header data ####
        data = sdf.read('./Data_n02/'+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        x  = data['Grid/Grid_mid'].data[0]/1.0e-6
        print('ok')
        y  = data['Grid/Grid_mid'].data[1]/1.0e-6
        X, Y = np.meshgrid(x, y) 
        name = '02'
        den_out = data['Derived/Number_Density/electron'].data/denunit
        den_out = den_out[250:dict_list[name],300-ywidth:300+ywidth]
        den_in = data['Derived/Number_Density/electronin'].data/denunit
        den_in = den_in[250:dict_list[name],300-ywidth:300+ywidth]
        ex = (data['Electric Field/Ey'].data-data['Electric Field/Ey_averaged'].data)/exunit
        ex = ex[250:dict_list[name],300-ywidth:300+ywidth]
        ex2 = ex**2

        plt.plot(y[300-ywidth:300+ywidth], den_out.sum(axis=0)/den_out[:,0].size, linewidth=3,linestyle='-',color='black', label=r'$n_{out}$')
        plt.plot(y[300-ywidth:300+ywidth], den_in.sum(axis=0)/den_in[:,0].size, linewidth=3,linestyle='-',color='red', label=r'$n_{in}$')
        plt.plot(y[300-ywidth:300+ywidth], 0.5*(den_in.sum(axis=0)/den_in[:,0].size+den_out.sum(axis=0)/den_out[:,0].size), linewidth=3,linestyle='--',color='orange', label=r'$\frac{n_{out}+n_{in}}{2}$')
        plt.plot(y[300-ywidth:300+ywidth], np.zeros_like(y[300-ywidth:300+ywidth])+0.2/2, linewidth=3,linestyle='--',color='green', label=r'$\frac{n_{ion}}{2}$')
        plt.xlabel('Y [$\lambda_0$]',fontdict=font)
        plt.ylabel(r'$n_e\ [n_c]$',fontdict=font, color='black')
        plt.ylim(0,1.1)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

        plt1 = plt.twinx()
        plt1.plot(y[300-ywidth:300+ywidth], (ex2.sum(axis=0)/ex2[:,0].size)**0.5 ,linewidth=3, linestyle='-',color='blue')
        plt1.set_ylabel(r'$\sqrt{\bar{E_y^2}}\ [m_ec\omega_0/e]$', color='blue')
        plt1.tick_params('y', colors='blue')
        plt1.set_ylim((0,15))
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);


        plt.subplot(2,3,4)
        #### header data ####
        data = sdf.read('./Data_n05/'+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        x  = data['Grid/Grid_mid'].data[0]/1.0e-6
        print('ok')
        y  = data['Grid/Grid_mid'].data[1]/1.0e-6
        X, Y = np.meshgrid(x, y) 
        name = '05'
        den_out = data['Derived/Number_Density/electron'].data/denunit
        den_out = den_out[250:dict_list[name],300-ywidth:300+ywidth]
        den_in = data['Derived/Number_Density/electronin'].data/denunit
        den_in = den_in[250:dict_list[name],300-ywidth:300+ywidth]
        ex = (data['Electric Field/Ey'].data-data['Electric Field/Ey_averaged'].data)/exunit
        ex = ex[250:dict_list[name],300-ywidth:300+ywidth]
        ex2 = ex**2

        plt.plot(y[300-ywidth:300+ywidth], den_out.sum(axis=0)/den_out[:,0].size, linewidth=3,linestyle='-',color='black', label=r'$n_{out}$')
        plt.plot(y[300-ywidth:300+ywidth], den_in.sum(axis=0)/den_in[:,0].size, linewidth=3,linestyle='-',color='red', label=r'$n_{in}$')
        plt.plot(y[300-ywidth:300+ywidth], 0.5*(den_in.sum(axis=0)/den_in[:,0].size+den_out.sum(axis=0)/den_out[:,0].size), linewidth=3,linestyle='--',color='orange', label=r'$\frac{n_{out}+n_{in}}{2}$')
        plt.plot(y[300-ywidth:300+ywidth], np.zeros_like(y[300-ywidth:300+ywidth])+0.5/2, linewidth=3,linestyle='--',color='green', label=r'$\frac{n_{ion}}{2}$')
        plt.xlabel('Y [$\lambda_0$]',fontdict=font)
        plt.ylabel(r'$n_e\ [n_c]$',fontdict=font, color='black')
        plt.ylim(0,1.1)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

        plt1 = plt.twinx()
        plt1.plot(y[300-ywidth:300+ywidth], (ex2.sum(axis=0)/ex2[:,0].size)**0.5 ,linewidth=3, linestyle='-',color='blue')
        plt1.set_ylabel(r'$\sqrt{\bar{E_y^2}}\ [m_ec\omega_0/e]$', color='blue')
        plt1.tick_params('y', colors='blue')
        plt1.set_ylim((0,15))
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

      
        plt.subplot(2,3,5)
        #### header data ####
        data = sdf.read('./Data_n10/'+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        x  = data['Grid/Grid_mid'].data[0]/1.0e-6
        print('ok')
        y  = data['Grid/Grid_mid'].data[1]/1.0e-6
        X, Y = np.meshgrid(x, y) 
        name = '10'
        den_out = data['Derived/Number_Density/electron'].data/denunit
        den_out = den_out[250:dict_list[name],300-ywidth:300+ywidth]
        den_in = data['Derived/Number_Density/electronin'].data/denunit
        den_in = den_in[250:dict_list[name],300-ywidth:300+ywidth]
        ex = (data['Electric Field/Ey'].data-data['Electric Field/Ey_averaged'].data)/exunit
        ex = ex[250:dict_list[name],300-ywidth:300+ywidth]
        ex2 = ex**2

        plt.plot(y[300-ywidth:300+ywidth], den_out.sum(axis=0)/den_out[:,0].size, linewidth=3,linestyle='-',color='black', label=r'$n_{out}$')
        plt.plot(y[300-ywidth:300+ywidth], den_in.sum(axis=0)/den_in[:,0].size, linewidth=3,linestyle='-',color='red', label=r'$n_{in}$')
        plt.plot(y[300-ywidth:300+ywidth], 0.5*(den_in.sum(axis=0)/den_in[:,0].size+den_out.sum(axis=0)/den_out[:,0].size), linewidth=3,linestyle='--',color='orange', label=r'$\frac{n_{out}+n_{in}}{2}$')
        plt.plot(y[300-ywidth:300+ywidth], np.zeros_like(y[300-ywidth:300+ywidth])+1.0/2, linewidth=3,linestyle='--',color='green', label=r'$\frac{n_{ion}}{2}$')
        plt.xlabel('Y [$\lambda_0$]',fontdict=font)
        plt.ylabel(r'$n_e\ [n_c]$',fontdict=font, color='black')
        plt.ylim(0,1.1)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

        plt1 = plt.twinx()
        plt1.plot(y[300-ywidth:300+ywidth], (ex2.sum(axis=0)/ex2[:,0].size)**0.5 ,linewidth=3, linestyle='-',color='blue')
        plt1.set_ylabel(r'$\sqrt{\bar{E_y^2}}\ [m_ec\omega_0/e]$', color='blue')
        plt1.tick_params('y', colors='blue')
        plt1.set_ylim((0,15))
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);


        plt.subplot(2,3,6)
        #### header data ####
        data = sdf.read('./Data_n20/'+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        x  = data['Grid/Grid_mid'].data[0]/1.0e-6
        print('ok')
        y  = data['Grid/Grid_mid'].data[1]/1.0e-6
        X, Y = np.meshgrid(x, y) 
        name = '20'
        den_out = data['Derived/Number_Density/electron'].data/denunit
        den_out = den_out[250:dict_list[name],300-ywidth:300+ywidth]
        den_in = data['Derived/Number_Density/electronin'].data/denunit
        den_in = den_in[250:dict_list[name],300-ywidth:300+ywidth]
        ex = (data['Electric Field/Ey'].data-data['Electric Field/Ey_averaged'].data)/exunit
        ex = ex[250:dict_list[name],300-ywidth:300+ywidth]
        ex2 = ex**2

        plt.plot(y[300-ywidth:300+ywidth], den_out.sum(axis=0)/den_out[:,0].size, linewidth=3,linestyle='-',color='black', label=r'$n_{out}$')
        plt.plot(y[300-ywidth:300+ywidth], den_in.sum(axis=0)/den_in[:,0].size, linewidth=3,linestyle='-',color='red', label=r'$n_{in}$')
        plt.plot(y[300-ywidth:300+ywidth], 0.5*(den_in.sum(axis=0)/den_in[:,0].size+den_out.sum(axis=0)/den_out[:,0].size), linewidth=3,linestyle='--',color='orange', label=r'$\frac{n_{out}+n_{in}}{2}$')
        plt.plot(y[300-ywidth:300+ywidth], np.zeros_like(y[300-ywidth:300+ywidth])+2.0/2, linewidth=3,linestyle='--',color='green', label=r'$\frac{n_{ion}}{2}$')
        plt.xlabel('Y [$\lambda_0$]',fontdict=font)
        plt.ylabel(r'$n_e\ [n_c]$',fontdict=font, color='black')
        plt.ylim(0,1.1)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);

        plt1 = plt.twinx()
        plt1.plot(y[300-ywidth:300+ywidth], (ex2.sum(axis=0)/ex2[:,0].size)**0.5 ,linewidth=3, linestyle='-',color='blue')
        plt1.set_ylabel(r'$\sqrt{\bar{E_y^2}}\ [m_ec\omega_0/e]$', color='blue')
        plt1.tick_params('y', colors='blue')
        plt1.set_ylim((0,15))
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);



        fig = plt.gcf()
        fig.set_size_inches(40, 13)
        fig.savefig('./cross_field_density'+str(n).zfill(4)+'.png',format='png',dpi=80)
        plt.close("all")


        print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

