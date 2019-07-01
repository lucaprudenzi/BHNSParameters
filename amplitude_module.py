import numpy as np

class Amplitude_class:

    def __init__(self, iterate=()):
        self.__dict__.update(iterate)
        self.m_1 = self.m_ns
        self.m_2 = self.m_bh
        self.chi_1 = self.chi_bh
        self.chi_2 = 0
    
    def Eta(self):
        eta = self.m_1*self.m_2/(self.m_1+self.m_2)**2
        return eta

    def Chi_s(self):
        chi_s = (self.chi_1*self.m_1+self.chi_2*self.m_2)/(self.m_1+self.m_2)
        return chi_s

    def Chi_a(self):
        chi_a = (self.chi_1*self.m_1-self.chi_2*self.m_2)/(self.m_1+self.m_2)
        return chi_a

    def Delta(self):
        delta = (1-4*self.Eta())**(1./2)
        return delta

    def A22220(self, real):
        a22220 = self.Eta()*(-0.6537*self.Chi_s()-4.0071) 
        return abs(a22220)

    def A21210(self, real):
        a21210 = self.Eta()*(2.3488*np.exp(2.6631j)*self.Delta()+\
                             0.8011*np.exp(5.7070j)*self.Chi_a()+\
                             3.5828*np.exp(5.5223j)*self.Eta()*self.Delta()+\
                             1.1744*np.exp(0.4254j)*self.Chi_s()*self.Delta()+\
                             0.6260*np.exp(5.3457j)*self.Chi_s()*self.Chi_a())
        if real:
            return a21210.real
        else:
            return abs(a21210)
    
    def A33330(self, real):
        a33330 = self.Eta()*(2.641*np.exp(2.9880j)*self.Delta()+\
                             1.6030*np.exp(0.6655j)*self.Delta()**2+\
                             1.0354*np.exp(3.6096j)*self.Chi_s()*self.Delta()+\
                             0.4911*np.exp(4.7347j)*self.Chi_a()**2)
        if real:
            return a33330.real
        else:
            return abs(a33330)

    def A32320(self, real):
        a32320 = self.Eta()*(2.5707*np.exp(4.1427j)*self.Eta()+\
                             9.4216*np.exp(0.8706j)*self.Eta()**2+\
                             0.5973*np.exp(2.1816j)*self.Eta()*self.Chi_s()+\
                             0.2104*np.exp(4.9043j)*self.Chi_a()**2+\
                             0.4417*np.exp(5.4543j)*self.Chi_a()*self.Delta()+\
                             0.9439*np.exp(1.7614)*self.Delta()**2)

        if real:
            return a32320.real
        else:
            return abs(a32320)

    def A32220(self, real):
        a32220 = self.Eta()*(1.3407*np.exp(2.9466j)*self.Eta()+\
                             0.0717*np.exp(5.5304j)+\
                             0.1061*np.exp(2.6432j)*self.Chi_s()**2+\
                             0.9894*np.exp(2.9294j)*self.Eta()*self.Chi_s()+\
                             0.3735*np.exp(3.3290j)*self.Chi_a()*self.Delta())
        if real:
            return a32220.real
        else:
            return abs(a32220)

    def A44440(self, real):
        a44440 = self.Eta()*(1.3284*np.exp(2.6831j)*self.Delta()**2+\
                             1.1619*np.exp(0.4142j)*self.Delta()**3+\
                             1.2790*np.exp(4.7226j)*self.Chi_s()*self.Chi_a()**2*self.Delta()+\
                             1.2387*np.exp(4.5616j)*self.Chi_s()*self.Chi_a()**3+\
                             1.2909*np.exp(2.8120j)*self.Chi_s()*self.Delta()**3+\
                             42.3575*np.exp(6.1418j)*self.Eta()**4)
        if real:
            return a44440.real
        else:
            return abs(a44440)
    
    def A43330(self,real):
        a43330 = self.Eta()*(0.0411*np.exp(2.6411j)*self.Chi_a()+\
                            0.0486*np.exp(3.2085j)*self.Chi_s()**2+\
                             0.8078*np.exp(2.7461j)*self.Eta()*self.Delta()+\
                             0.1940*np.exp(3.0292j)*self.Chi_s()*self.Delta()+\
                             0.0529*np.exp(3.5830j)*self.Chi_a()**2+\
                             0.0358*np.exp(0.1731j)*self.Delta()**2)
        if real:
            return a43330.real
        else:
            return abs(a43330)

    def A43430(self,real):
        a43430 = self.Eta()*(0.5665*np.exp(3.3993j)*self.Delta()+\
                             0.1457*np.exp(4.7476j)*self.Chi_a()+\
                             0.8239*np.exp(1.8174j)*self.Eta()*self.Chi_a()+\
                             0.0507*np.exp(4.7495j)*self.Chi_a()*self.Chi_s()+\
                             0.9806*np.exp(0.6029j)*self.Delta()**3+\
                             10.1678*np.exp(6.2185j)*self.Eta()**2*self.Delta())
        if real:
            return a43430.real
        else:
            return abs(a43430)

