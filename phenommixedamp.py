import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import fsolve
import spin_module as sm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import lal

class PhenomBHNS_class:

    xi01=[-0.0024172093701738235,0.5962010933576927,-32.83285391505444,161.89021265549897,-631.9871537581226,-48.086652464073154,4.148542851248616,-0.05471995193663809,-1.2349578240136077]
    xi02=[-0.0010926532157647696,-0.05600496937966142,8.859039052577113,-47.0230352671574,246.34764324446115,-364.31302921602577,-4.069743432687181,0.020944687428888942,0.34233508992115164]
    xi11=[-0.019168568675124984,0.15195931833090354,29.30762955550024,-175.05034482339252,1047.846454700688,-521.4742751778888,-87.52069774952481,0.3553835917253634,6.0621729120758765]
    xi10=[0.07267102554435025,-2.9698758936751806,79.54388527811831,-322.46920438601825,335.48879536851206,1869.6810123048372,-48.97027656287295,0.11513830998008491,5.949166424824434]
    xi20=[-0.25037513499482245,13.119790089123903,-434.8761203428741,1586.8887503860365,-5114.768883285114,735.377346283158,666.4552087394841,0.9639539830624523,-10.693356457030575]
    
    def __init__(self, iterate=()):
        self.__dict__.update(iterate)
        
        # useful quantities
        self.eta = self.q/(1+self.q)**2
        self.m_bh = self.q*self.m_ns
        self.M = self.m_ns+self.m_bh
        #self.G = 6.67*10**(-8)
        self.G = lal.G_SI
        #self.c = 3*10**(10)
        self.c = lal.C_SI
        #self.m_sol = 2*10**33
        self.m_sol = lal.MSUN_SI
        self.chi_ns=0
        self.winfact = 1./2 # ok=1/2, paper=1
        # coefficients for PhenomC model
        self.gamma = self.xi01[6]*self.chi_bh+self.xi02[6]*self.chi_bh**2+self.xi11[6]*self.eta*self.chi_bh+self.xi10[6]*self.eta+self.xi20[6]*self.eta**2
        self.delta1 = self.xi01[7]*self.chi_bh+self.xi02[7]*self.chi_bh**2+self.xi11[7]*self.eta*self.chi_bh+self.xi10[7]*self.eta+self.xi20[7]*self.eta**2
        self.delta2 = self.xi01[8]*self.chi_bh+self.xi02[8]*self.chi_bh**2+self.xi11[8]*self.eta*self.chi_bh+self.xi10[8]*self.eta+self.xi20[8]*self.eta**2
        # spin section
        binary_par = {
            "m_ns": self.m_ns,
            "q": self.q,              
            "theta_i_deg": self.theta_i_deg, 
            "chi_bh": self.chi_bh,      
            "EOS": self.EOS,         
            "type": self.type
        }

        spin_object = sm.Spin_class(binary_par)
        chi_f, theta_f, m_f, m_out = spin_object.FinalPar()
        self.chi_f = float(chi_f)
        self.theta_f = float(theta_f)
        self.m_f = float(m_f)
        self.m_out = float(m_out)
        
        # Ringdown frequency and damping time
        # From Arxiv:gr-qc/0512160v2 eq. E1, E2
        self.f_220_ringdown = (1.5251-1.1568*(1-self.chi_f)**(0.1292))/(2*np.pi*self.m_f*self.m_sol*self.G/self.c**2)*self.c
        self.Q_220_ringdown = 0.7+1.4187*(1-self.chi_f)**(-0.4990)
        self.tau_220_ringdown = self.Q_220_ringdown/self.f_220_ringdown/np.pi

        # PhenomC for BBH
        if self.type == "BBH":
            self.eps_insp = 1
            self.k = 1
            self.sigma_tid = 0
            self.eps_tid = 1
            
            self.f0_pn = 0.98*self.f_220_ringdown
            self.f0_pm = 0.98*self.f_220_ringdown
            self.f0_rd = 0.98*self.f_220_ringdown
        
        # PhenomMixed
        if self.type == "BHNS":
            print("Total mass outside BH at late times: %.3f" %m_out)
            m_g, m_b, compactness = np.loadtxt("EOS/equil_{}.d".format(self.EOS), unpack=True)
            C_from_mg = interp1d(m_g, compactness)
            self.C = float(C_from_mg(self.m_ns))
            mb_from_mg = interp1d(m_g, m_b)
            self.m_ns_b = float(mb_from_mg(self.m_ns))
            
            xi_tid = fsolve(spin_object.XiTidal,10)
            r_tid = xi_tid*(self.m_ns*self.m_sol*self.G/self.c**2)/self.C*(1-2*self.C)
            #r_tid = self.eta**(-1./3)*(self.m_ns*self.m_sol*self.G/self.c**2)/self.C*(1-2*self.C)

            # frequencies in Hz
            self.f_tid = float(1./(np.pi*(self.chi_f*self.m_f*self.m_sol*self.G/self.c**2+(r_tid**3/(self.m_f*self.m_sol*self.G/self.c**2))**(1./2)))*self.c)
            #print("theta_i_deg", self.theta_i_deg)
            #print("chi_bh", self.chi_bh)
            #print("chi_f", self.chi_f)
            self.f_220_tilde = 0.99*0.98*self.f_220_ringdown
            self.k = 1.25
            if self.x_var == "f":
                plt.axvline(self.f_220_ringdown,0.9,1, color="c")
                plt.axvline(self.f_tid, 0.9, 1, color="r")
            if self.x_var == "fm":
                plt.axvline(self.f_220_ringdown*self.G/self.c**3*self.M*self.m_sol,0.9,1, color="c")
                plt.axvline(self.f_tid*self.G/self.c**3*self.M*self.m_sol, 0.9, 1, color="r")
            if self.eta >= 5./36:            
                self.alpha = 1+self.winfact*self.Windowed_function(self.eta, 0.146872, 0.0749456, -1)
            if self.eta < 5./36:
                self.alpha = 1+self.winfact*self.Windowed_function(5./36, 0.146872, 0.0749456, -1)
            self.delta2prime = self.alpha*self.delta2

            # Non disruptive
            if self.f_220_ringdown <= self.f_tid and self.m_out == 0:
                print("Case: Non disruptive")
                self.x_nd1 = ((self.f_tid-self.f_220_tilde)/self.f_220_tilde)**2-0.571505*self.C-0.00508451*self.chi_bh
                self.x_nd2 = ((self.f_tid-self.f_220_tilde)/self.f_220_tilde)**2-0.657424*self.C-0.0259977*self.chi_bh
                
                self.eps_tid = self.winfact*2*self.Windowed_function(self.x_nd1,-0.0796251,0.0801192,1)
                self.sigma_tid = self.winfact*2*self.Windowed_function(self.x_nd2,-0.206465,0.226844,-1)/(self.G/self.c**3*self.M*self.m_sol)
                #self.delta2prime = 1.62496*self.winfact*self.Windowed_function((self.f_tid-self.f_220_tilde)/self.f_220_tilde, 0.0188092, 0.338737, -1)
                
                self.f0_pn = self.f_220_tilde
                self.f0_pm = self.f_220_tilde
                self.f0_rd = self.f_220_tilde
                #print("")
                #print("f0_pn", self.f0_pn)
                #print("f0_pm", self.f0_pm)
                #print("f0_rd", self.f0_rd)
                #print("eps_tid", self.eps_tid)
                #print("sigma_tid", self.sigma_tid)
                #print("delta2prime", self.delta2prime)
                #print("x_nd", self.x_nd1)
                #print("")


            # Disruptive
            if self.f_tid < self.f_220_ringdown and self.m_out > 0:
                print("Case: Disruptive")
                self.x_d1 = self.m_out/self.m_ns_b+0.424912*self.C+\
                        0.363604*(self.eta)**(1./2)-0.0605591*self.chi_bh
                self.x_d2 = self.m_out/self.m_ns_b-0.132754*self.C+\
                        0.576669*(self.eta)**(1./2)-0.0603749*self.chi_bh-\
                        0.0601185*self.chi_bh**2-0.0729134*self.chi_bh**3
                
                self.esp_tid = 0
                self.eps_insp = -1.61724*self.x_d1+1.29971
                self.sigma_tid = (-0.293237*self.x_d2+0.1377221)/\
                        (self.G/self.c**3*self.M*self.m_sol)
                
                self.f0_pn = self.eps_insp*self.f_tid
                self.f0_pm = self.f_tid
            
            # Mildly disruptive 1
            if self.f_tid < self.f_220_ringdown and self.m_out == 0:
                print("Case: Midly disruptive, no ringdown, no M outside BH")
                self.x_d1 = self.m_out/self.m_ns_b+0.424912*self.C+\
                        0.363604*(self.eta)**(1./2)-0.0605591*self.chi_bh
                self.x_d2 = self.m_out/self.m_ns_b-0.132754*self.C+\
                        0.576669*(self.eta)**(1./2)-0.0603749*self.chi_bh-\
                        0.0601185*self.chi_bh**2-0.0729134*self.chi_bh**3
                self.x_nd2 = ((self.f_tid-self.f_220_tilde)/self.f_220_tilde)**2-0.657424*self.C-0.0259977*self.chi_bh
                
                self.eps_tid = 0
                self.eps_insp = -1.61724*self.x_d1+1.29971
                self.sigma_tid_d = (-0.293237*self.x_d2+0.1377221)/\
                        (self.G/self.c**3*self.M*self.m_sol)
                self.sigma_tid_nd = self.winfact*2*self.Windowed_function(self.x_nd2,-0.206465,0.226844,-1)/(self.G/self.c**3*self.M*self.m_sol)
                self.sigma_tid = (self.sigma_tid_d+self.sigma_tid_nd)/2.
                
                self.f0_pn = (1-self.q**(-1))*self.f_220_tilde+\
                        self.q**(-1)*self.eps_insp*self.f_tid
                self.f0_pm = (1-self.q**(-1))*self.f_220_tilde+\
                        self.q**(-1)*self.f_tid

            # Mildly disruptive 2
            if self.f_220_ringdown <= self.f_tid and self.m_out > 0:
                print("Case: Mildly diruptive, ringdown, M outside BH")
                self.x_d1 = self.m_out/self.m_ns_b+0.424912*self.C+\
                        0.363604*(self.eta)**(1./2)-0.0605591*self.chi_bh
                self.x_nd1 = ((self.f_tid-self.f_220_tilde)/self.f_220_tilde)**2-0.571505*self.C-0.00508451*self.chi_bh
                self.x_nd2 = ((self.f_tid-self.f_220_tilde)/self.f_220_tilde)**2-0.657424*self.C-0.0259977*self.chi_bh
                
                self.eps_insp = -1.61724*self.x_d1+1.29971
                self.eps_tid = self.winfact*2*self.Windowed_function(self.x_nd1,-0.0796251,0.0801192,1)
                self.sigma_tid = self.winfact*2*self.Windowed_function(self.x_nd2,-0.206465,0.226844,-1)/(self.G/self.c**3*self.M*self.m_sol)
                #self.delta2prime = 1.62496*self.winfact*self.Windowed_function((self.f_tid-self.f_220_tilde)/self.f_220_tilde, 0.0188092, 0.338737, -1)
                self.f0_pn = self.eps_insp*self.f_220_tilde
                self.f0_pm = self.f_220_tilde
                self.f0_rd = self.f_220_tilde

            print("Tidal frequency %.2f" %self.f_tid)
            print("Ringdown frequency: %.2f" %self.f_220_tilde)
        #print("eps_tid", self.eps_tid)
        #print("eps_insp", self.eps_insp)
        #print("m_out", self.m_out)
        #print("f_tidal", self.f_tid)
        #print("f_220", self.f_220_tilde)
        #print("f_220_ring", self.f_220_ringdown)

        self.d_bbh = 0.015/(self.G/self.c**3*self.M*self.m_sol)
        self.d = self.d_bbh+self.sigma_tid

    def Windowed_function(self, f, f0, d, sign):
        w = 1./2*(1+sign*np.tanh(4*(f-f0)/d))
        return w

    def Lorentzian(self, f, f0, sigma):
        lor = 1./(1./4+((f-f0)/sigma)**2)
        return lor


    # Amplitude in time domain (manca D_l al denominatore ma si plotta A*D_l) 
    def Amp_PN_T(self,f):
        x = (np.pi*f*self.G/self.c**3*self.M*self.m_sol)**(2./3)
        chi_bh = self.chi_bh
        chi_ns = self.chi_ns
        chi=(chi_bh*self.m_bh+chi_ns*self.m_ns)/self.M
        eta = self.eta

        amp_t = 8*eta*x*(np.pi/5.)**(1./2)*\
                 (1+\
                  x*(-107./42.+55*eta/42.)+\
                  x**(3./2.)*(2*np.pi-4./3*chi+2*eta/3.*(chi_bh+chi_ns))+\
                  x**2*(-2173./1512.-1069*eta/216.+2047*eta**2/\
                        1512.+2*eta*chi_bh*chi_ns)+\
                  x**(5./2.)*(-107*np.pi/21.-24.*1j*eta+\
                              34*np.pi*eta/21.)+\
                  x**3*(27027409./646800.-856*np.euler_gamma/105.+\
                         428*1j/105.*np.pi+2*np.pi**2/3.+\
                         (-278185./33264.+41*np.pi**2/96.)*\
                         eta-20261*eta**2/2772.+\
                         114635*eta**3/99792.-428*np.log(16*x)/105.))
        
        return amp_t

    def TaylorT4xdot(self,f):
        
        x = (np.pi*f*self.G/self.c**3*self.M*self.m_sol)**(2./3)
        chi = (self.chi_bh*self.m_bh+self.chi_ns*self.m_ns)/self.M
        chi_ns = self.chi_ns
        chi_bh = self.chi_bh
        eta = self.eta

        xdot = 64./5.*self.eta*x**5*\
                (1+\
                 x*(-743./336-11./4.*self.eta)+\
                 x**(3./2)*(4*np.pi-113./12*chi+19./6.*eta*chi_ns+\
                            chi_bh)+\
                 x**2*(34103./18144+5*chi**2+eta*\
                       (13661./2016-chi_bh*chi_ns/8.)+\
                       59./18.*eta**2)+\
                 x**(5./2)*(-np.pi*(4159./672.+189./8.*eta)-\
                            chi*(31571./1008-1165./24.*eta)+\
                            (chi_ns+chi_bh)*(21863./1008.*eta-76./6.*eta**2)-\
                            3./4.*chi**3+9./4.*eta*chi*chi_bh*chi_ns)+\
                 x**3*(16447322263./139708800-1712./105*np.euler_gamma+\
                       16./3.*np.pi**2-856./105*np.log(16*x)+\
                       eta*(451./48.*np.pi**2-56198689./217728.)+\
                       541./896*eta**2-5605./2592.*eta**3-\
                       80./3*np.pi*chi+(20./3*np.pi-1135./36*chi)*eta*\
                       (chi_bh+chi_ns)+(64153./1008-457./36*eta)*chi**2-\
                       (787./144*eta-3037./144*eta**2)*chi_bh*chi_ns)+
                 x**(7/2)*(-np.pi*(4415./4032.-358675./6048*eta-\
                                   91495./1512*eta**2)-\
                           chi*(2529407./27216-845827./6048*eta+\
                                41551./864*eta**2)+\
                           (chi_bh+chi_ns)*(1580239./54432*eta-451597./6048*\
                                            eta**2+2045./432.*eta**3+107./6.*\
                                            eta*chi**2-5./24.*eta**2*\
                                            chi_bh*chi_ns)+\
                           12*np.pi*chi**2-chi**3*(1505./24.+1./8.*eta)+\
                           chi*chi_bh*chi_ns*(101./24.*eta+3./8*eta**2)))
        return xdot

    # amplitude in frequency domain 
    def Amp_PN_F(self,f):
        prefactor = self.G/self.c**3*self.M*self.m_sol
        factor_220 =  (5/(4*np.pi))**(1./2)/2.
        factor_220 = 1
        amp_f = self.Amp_PN_T(f)*(2.*np.pi/(3*(np.pi*f)**(1./3.)*prefactor**(1./3)*self.TaylorT4xdot(f)))**(1./2)*factor_220
        
        return amp_f

    def Amp_PM(self, f):
        prefactor = self.G/self.c**3*self.M*self.m_sol
        a_pm = self.k*self.gamma*(f*prefactor)**(5./6)
        return a_pm
    
    def Amp_RD(self, f):
        prefactor = self.G/self.c**3*self.M*self.m_sol
        if self.type == "BBH":
            lor = self.Lorentzian(f, self.f0_rd, self.delta2/(np.pi*self.tau_220_ringdown))
        if self.type== "BHNS":
            lor = self.Lorentzian(f, self.f0_rd, self.delta2prime/(np.pi*self.tau_220_ringdown))
        a_rd = self.eps_tid*self.delta1*lor*(f*prefactor)**(-7./6)
        return a_rd

    def Amp_Phenom(self,f):
        a_pn = self.Amp_PN_F(f)
        w_pn = self.Windowed_function(f, self.f0_pn, self.d, -1)

        a_pm = self.Amp_PM(f)
        w_pm = self.Windowed_function(f, self.f0_pm, self.d, -1)
        factor_220 = (5./(np.pi*4))**(1./2)/2.

        if self.type == "BBH":
            a_rd = self.Amp_RD(f)
            w_rd = self.Windowed_function(f, self.f0_rd, self.d, 1)
            a_phenom = a_pn*w_pn+a_pm*w_pm+a_rd*w_rd
            #a_phenom = a_pm*w_pm+a_rd*w_rd

        if self.type == "BHNS":
            if self.f_tid < self.f_220_ringdown:
                a_rd = 0
                w_rd = 0
            if self.f_220_ringdown <= self.f_tid:
                a_rd = self.Amp_RD(f)
                w_rd = self.Windowed_function(f, self.f0_rd, self.d, 1)
            a_phenom = a_pn*w_pn+a_pm*w_pm+a_rd*w_rd

        return a_phenom*factor_220
