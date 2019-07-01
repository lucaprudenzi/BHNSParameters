import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
 
class Detector_class:
    def __init__(self, iterate=()):
        self.__dict__.update(iterate)

    def SNR(self, f_arr, amp, detector):
        if detector == "LIGO":
            psd = self.PSD_LIGO(f_arr)
        if detector == "Virgo":
            psd = self.PSD_Virgo(f_arr)
        if detector == "ET":
            psd = self.PSD_ET(f_arr)
        df = f_arr[2]-f_arr[1]
        snr = 0   
        #f_arr2 = [i for i in f_arr if i<40]

        for j in range(len(f_arr)):
            snr = snr+(amp[j]/f_arr[j])**2/psd[j]**2*df
        snr = (4*snr)**(1./2)
        return snr 
    
    def Plot(self, f_arr, detector):
        if detector == "LIGO":
            psd = self.PSD_LIGO(f_arr)
        if detector == "Virgo":
            psd = self.PSD_Virgo(f_arr)
        if detector == "ET":
            psd = self.PSD_ET(f_arr)
#       psd1 = self.PSD_Virgo(f_arr)
        plt.plot(f_arr, psd, label=r"$\sqrt{S_h(f)} %s$" %detector)
#       plt.plot(f_arr,psd1,label="Virgo")

    def PSD_LIGO(self, f_arr):
        f, psd = np.loadtxt("PSD/LIGO-P1200087-v18-aLIGO_DESIGN.txt", unpack=True)
        asd_from_f = interp1d(f, psd)
        psd_arr=[]
        for i in range(len(f_arr)):
            psd_arr.append(asd_from_f(f_arr[i]))

        #f0 = 215. # Hz
        #x = f/f0
        #psd_arr = 10**(-49)*(x**(-4.14)-5*x**(-2)+111*(1-x**2+x**4/2.)/(1+x**2/2.))
        return psd_arr

    def PSD_Virgo(self, f_arr):
        f, psd = np.loadtxt("PSD/LIGO-P1200087-v18-AdV_DESIGN.txt", unpack=True)
        psd_from_f = interp1d(f, psd)
        psd_arr=[]
        for i in range(len(f_arr)):
            psd_arr.append(psd_from_f(f_arr[i]))

        #f0 = 720. # Hz
        #x = f/f0
        #psd_arr = 10**(-47)*(2.67*10**(-7)*x**(-5.6)+0.59*np.exp((np.log(x))**2*(-3.2-1.08*np.log(x)-0.13*(np.log(x))**2))*x**(-4.1)+0.68*np.exp(-0.73*(np.log(x))**2)*x**5.34)
        return psd_arr
    
    def PSD_ET(self, f):
        f0 = 100.
        x = f/f0
        psd = 2.39*10**(-27)*x**(-15.64)+0.349*x**(-2.145)+1.76*x**(-0.12)+0.409*x**(1.10)
        return psd
