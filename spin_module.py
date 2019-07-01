import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math 

class Spin_class:

    # Constants
    deg_rad = np.pi/180.
    rad_deg = 180./np.pi

    def __init__(self, iterate=()):
        self.__dict__.update(iterate)
        self.theta_i = self.theta_i_deg*self.deg_rad
        self.m_bh = self.q*self.m_ns
        self.M = self.m_bh+self.m_ns
        self.nu = self.m_bh*self.m_ns/(self.m_bh+self.m_ns)**2 # symmetric mass ratio
        self.eta = self.q/(1+self.q)**2
        m_g, m_b, compactness = np.loadtxt("EOS/equil_{}.d".format(self.EOS), unpack=True)
        C_from_mg = interp1d(m_g, compactness)
        self.C = C_from_mg(self.m_ns)
        mb_from_mg = interp1d(m_g, m_b)
        self.m_ns_b = mb_from_mg(self.m_ns)
        self.r_ns = self.m_ns/self.C # M -> G/c^2 M
        
        #m_ns_b_test = self.M_ns_bar()
        #print("m_test", m_ns_b_test)

    # ISCO radius in mass unit
    def R_isco(self,chi,direction):
        #print("Risco is using a", a)
        #print(2+(1-a**2)**(1./3)*((1+a)**(1./3)+(1-a)**(1./3)))
        #print("")
        if math.isnan(chi):
            print(chi) 
        Z1 = 1+(1-chi**2)**(1./3)*((1+chi)**(1./3)+(1-chi)**(1./3))
        Z2 = (3*chi**2+Z1**2)**(1./2)
        r_isco = (3+Z2-direction*((3-Z1)*(3+Z1+2*Z2))**(1./2))
        return r_isco
    
    # Energy for test particle for equatorial orbit
    # Spin and angular momentum are along the same direction
    def E_aligned(self, r, chi, direction):
        e = (r**2-2*r+direction*chi*r**(1./2))/(r*(r**2-3*r+direction*2*chi*r**(1./2))**(1./2))
        return e

    def E(self,chi,theta):
        r_isco_pro = self.R_isco(np.abs(chi), 1)
        r_isco_ret = self.R_isco(np.abs(chi), -1)
        e_angle_pro = 1./2*(1+np.cos(theta))*self.E_aligned(r_isco_pro,np.abs(chi),1)
        e_angle_ret = np.abs(1./2*(1-np.cos(theta))*self.E_aligned(r_isco_ret,np.abs(chi),-1))
        e_angle = e_angle_pro+e_angle_ret
        return e_angle
        
    def E_rad(self):
        e_i = self.E(self.chi_bh, self.theta_i)
        e_rad = self.M*(1-e_i)*self.nu
        return e_rad

    # orbital angular momentum per unit mass of a test particle 
    # [equatorial orbit]
    # spin and angular momentum aligned (along z direction)
    def L_z(self, r, chi, direction):
        l_z = direction*(r**2-direction*2*chi*r**(1./2)+chi**2)/(r**(1./2)*(r**2-3*r+direction*2*chi*r**(1./2))**(1./2))
        return l_z

    # Spin and angular momentum not aligned
    def L_orb(self,chi,theta):
        
        r_isco_pro = self.R_isco(np.abs(chi), 1)
        r_isco_ret = self.R_isco(np.abs(chi), -1)
        l_orb_pro = 1./2*(1+np.cos(theta))*self.L_z(r_isco_pro,np.abs(chi),1)
        l_orb_ret = np.abs(1./2*(1-np.cos(theta))*self.L_z(r_isco_ret,np.abs(chi),-1))
        l_orb = l_orb_pro+l_orb_ret
        return l_orb

    # baryon mass of the neutron star
    def M_ns_bar(self):
        baryonic_mass = self.m_ns*(1+0.6*self.C/(1-0.5*self.C)) 
        return baryonic_mass

    def R_isco_eff(self):
        theta_i = self.theta_i
        chi_bh = self.chi_bh
        r_isco_angle = 1./2*(1+np.cos(theta_i))*self.R_isco(np.abs(chi_bh),1)+1./2*(1-np.cos(theta_i))*self.R_isco(np.abs(chi_bh),-1)
        return r_isco_angle

    def Chi_eff(self, chi_eff):
        theta_i = self.theta_i
        chi_bh = self.chi_bh
        if chi_eff>1:
            chi_eff = 0.999
        if chi_eff<0:
            chi_eff = 0
        chi_eff = np.abs(chi_eff)
        r_isco_angle = 1./2*(1+np.cos(theta_i))*self.R_isco(np.abs(chi_bh),1)+1./2*(1-np.cos(theta_i))*self.R_isco(np.abs(chi_bh),-1)
        #print("r_isco_angle", r_isco_angle)
        r_isco_noangle = self.R_isco(np.abs(chi_eff),1)
        #print("r_isco_noangle", r_isco_noangle)
        f = r_isco_angle-r_isco_noangle
        
        return f

    def XiTidal(self,xi):
        if math.isnan(xi):
            print(xi)
        m_bh = self.m_bh
        m_ns = self.m_ns
        chi_bh = self.chi_bh
        k = self.q*self.C
        
        f = m_ns/m_bh*xi**3-3*(xi**2-2*k*xi+chi_bh**2*k**2)/(xi**2-3*k*xi+2*chi_bh*(k**3*xi)**(1./2))
        return f
   
    def M_out2012(self):
        if self.type == "BBH":
            m_out = 0
        if self.type == "BHNS":
            # Foucart 2012 relativistic
            chi_eff = fsolve(self.Chi_eff,0.5) 
            xi_tidal = fsolve(self.XiTidal,10)
            #print("q", self.q)
            #print("chi_bh", self.chi_bh)
            #print("chi_eff", chi_eff)
            #print("xi",xi_tidal)
            #print("")
            r_isco = self.R_isco(np.abs(chi_eff),1)
            alpha = 0.296
            beta = 0.171
            m_out = self.m_ns_b*(alpha*xi_tidal*(1-2*self.C)-beta*r_isco*self.q*self.C)
            if m_out<0:
                m_out=0
            
            # TEST
            #m_out = self.m_ns_b*(0.288*(3*self.q)**(1./3)*(1-2*self.C)-0.148*r_isco*self.q*self.C)
            #d_tidal2018 =  self.r_ns*self.eta**(-1./3)*(1-2*self.C)
            #d_tidal2012 = self.r_ns*xi_tidal*(1-2*self.C)
            #f2012 = 1./(d_tidal2012**3/(self.m_ns+self.m_bh))**(1./2)
            #f2018 = 1./(d_tidal2018**3/(self.m_ns+self.m_bh))**(1./2)
            #delta = f2018-f2012
            #print("q",self.q)
            #print("chi_bh", self.chi_bh)
            #print("delta", delta)
            #print("")
              
        return m_out

    # mass outside the black hole at late times
    def M_out(self):
        if self.type == "BBH":
            m_out = 0
        if self.type == "BHNS":
            m_bh = self.m_bh
            q = self.q
            C = self.C
            eta = self.eta
 
            # Foucart 2018
            r_isco_eff = self.R_isco_eff()
            alpha = 0.406
            beta = 0.139
            gamma = 0.255
            delta = 1.761

            m_out = (alpha*(1-2*C)*eta**(-1./3)-beta*r_isco_eff*C*eta**(-1)+gamma)
            if m_out<0:
                m_out = 0
            else:
                m_out = m_out**(delta)
            m_out = m_out*self.m_ns_b
        
        return m_out

    # mass fraction that does tidal distruption
    def Fraction(self):
        fraction = 0
        nu = self.nu

        if (nu<=0.16):
            fraction = 0
        if (nu>0.16 and nu<2./9):
            fraction = 1./2*(1-np.cos(np.pi*(nu-0.16)/(2./9-0.16)))
        if (nu>2./9):
            fraction = 1
        
        return fraction
        
    # Mass of the final black hole
    def M_f(self, chi_f, theta_f):
        e_rad = self.E_rad()
        e_f = self.E(chi_f, theta_f)
        m_out = self.m_out
        m_f = self.M-e_rad-e_f*m_out

        return m_f
   
    # function to find a_final for third Pannarale formula (generic spin)
    def ChifThetaf(self,p):
        chi_f, theta_f = p
        if chi_f>=1:
            chi_f = 0.99999999
        if chi_f<=-1:
            chi_f = -0.99999
        
        if self.type == "BHNS":
            chi_f = np.abs(chi_f) 
            l_orb = self.L_orb(chi_f, theta_f)
            m_f = self.M_f(chi_f,theta_f)
            m_out = self.m_out
            fraction = self.Fraction()
            m_ns_b = self.m_ns_b
        
            m_ns = self.m_ns
            m_bh = self.m_bh
            chi_bh = self.chi_bh
            theta_i = self.theta_i
            
            f1 = chi_f*m_f**2-(chi_bh*m_bh**2*np.cos(theta_i-theta_f)+l_orb*m_bh*((1-fraction)*m_ns+fraction*m_ns_b-m_out)*np.cos(theta_f))

            f2 = chi_bh*m_bh**2*np.sin(theta_i-theta_f)-l_orb*m_bh*((1-fraction)*m_ns+fraction*m_ns_b-m_out)*np.sin(theta_f)
        
        if self.type == "BBH":
            chi_f = np.abs(chi_f) 
        
            l_orb = self.L_orb(chi_f, theta_f)
            m_f = self.M_f(chi_f,theta_f)
            
            m_bh2 = self.m_ns
            m_bh1 = self.m_bh
            chi_bh1 = self.chi_bh
            theta_i = self.theta_i

            f1 = chi_f*m_f**2-(chi_bh1*m_bh1**2*np.cos(theta_i-theta_f)+l_orb*m_bh1*m_bh2*np.cos(theta_f))

            f2 = chi_bh1*m_bh1**2*np.sin(theta_i-theta_f)-l_orb*m_bh1*m_bh2*np.sin(theta_f)

        return (f1,f2)
    
    def ChifAligned(self,chi_f):
        if chi_f>=1:
            chi_f = 0.99999999
        if chi_f<=-1:
            chi_f = -0.999999
        if self.theta_i_deg == 0:
            sign = 1
        if self.theta_i_deg == 180:
            sign = -1

        if self.type == "BHNS":
            r_isco_f = self.R_isco(chi_f, sign)
            l_orb = self.L_z(r_isco_f, chi_f, sign)
            m_out = self.m_out
            fraction = self.Fraction()
            m_ns_b = self.m_ns_b 

            e_rad = self.E_rad()
            e_f = self.E_aligned(r_isco_f, chi_f, sign)
            m_out = self.m_out
            m_f = self.M-e_rad-e_f*m_out
       
            m_ns = self.m_ns
            m_bh = self.m_bh
            chi_bh = self.chi_bh
            theta_i = self.theta_i
            
            f = chi_f*m_f**2-(chi_bh*m_bh**2+l_orb*m_bh*((1-fraction)*m_ns+fraction*m_ns_b-m_out))

        if self.type == "BBH":
            r_isco_f = self.R_isco(chi_f, sign)
            l_orb = self.L_z(r_isco_f, chi_f, sign)
            
            e_rad = self.E_rad()
            e_f = self.E_aligned(r_isco_f, chi_f, sign)
            m_out = self.m_out
            m_f = self.M-e_rad-e_f*m_out

            m_bh2 = self.m_ns
            m_bh1 = self.m_bh
            chi_bh1 = self.chi_bh
            theta_i = self.theta_i

            f = chi_f*m_f**2-(chi_bh1*m_bh1**2+l_orb*m_bh1*m_bh2)
        return f

    def FinalPar(self):
        self.m_out = float(self.M_out())
        #self.m_out = float(self.M_out2012())


        # initial conditions of the final parameters
        chi_f_guess = 0
        theta_f_guess = 0
        
        # fine tuning of theta_f_guess initial condition 
        if self.theta_i_deg == 0:
            chi_f = fsolve(self.ChifAligned, chi_f_guess, maxfev=10000)
            theta_f = 0
            m_f = self.M_f(chi_f,theta_f)
            return chi_f, theta_f, m_f, self.m_out

        if self.theta_i_deg > 0 and self.theta_i_deg <= 140:
            theta_f_guess = np.pi/4
        if self.theta_i_deg > 140 and self.theta_i_deg<=172:
            if self.chi_bh<0.4:
                theta_f_guess=np.pi/4
            if self.chi_bh>=0.4:
                theta_f_guess = np.pi/2-0.10# in rad
        if self.theta_i_deg>172 and self.theta_i<=178:
            if self.chi_bh<0.5:
                if self.q<=5:
                    theta_f_guess = 45*np.pi/180. 
                if self.q>5:
                    theta_f_guess = 85*np.pi/180.
            if self.chi_bh>=0.5 and self.chi_bh<=0.6:
                if self.q<=5:
                    theta_f_guess = 45*np.pi/180. 
                if self.q>5:
                    theta_f_guess = 92*np.pi/180.
            if self.chi_bh>0.6:
                if self.q<=5:
                    theta_f_guess = 90*np.pi/180
                if self.q>5:
                    theta_f_guess = 150*np.pi/180.
        
        if self.theta_i_deg>178 and self.theta_i_deg<=179:
            if self.chi_bh<0.4:
                if self.q<=5:
                    theta_f_guess = 45*np.pi/180. 
                if self.q>5:
                    theta_f_guess = 85*np.pi/180.
            if self.chi_bh>=0.4 and self.chi_bh<0.5:
                if self.q<=5:
                    theta_f_guess = 45*np.pi/180. 
                if self.q>5:
                    theta_f_guess = 87*np.pi/180.
            if self.chi_bh>=0.5 and self.chi_bh<=0.6:
                if self.q<=5:
                    theta_f_guess = 45*np.pi/180. 
                if self.q>5:
                    theta_f_guess = 93*np.pi/180.
            if self.chi_bh>0.6 and self.chi_bh<0.9:
                if self.q<=5:
                    theta_f_guess = 86*np.pi/180
                if self.q>5:
                    theta_f_guess = 150*np.pi/180.
            if self.chi_bh>=0.9:
                if self.q<=4:
                    theta_f_guess = 85*np.pi/180
                if self.q>=4 and self.q<6:
                    theta_f_guess = 150*np.pi/180.
        
        if self.theta_i_deg>179 and self.theta_i_deg<179.3:
            if self.chi_bh<0.4:
                if self.q<=5:
                    theta_f_guess = 45*np.pi/180. 
                if self.q>5:
                    theta_f_guess = 85*np.pi/180.
            if self.chi_bh>=0.4 and self.chi_bh<0.5:
                if self.q<=7:
                    theta_f_guess = 20*np.pi/180. 
                if self.q>7:
                    theta_f_guess = 90*np.pi/180.
            if self.chi_bh>=0.5 and self.chi_bh<0.6:
                if self.q<=5:
                    theta_f_guess = 30*np.pi/180. 
                if self.q>5 and self.q<=7:
                    theta_f_guess = 80*np.pi/180.
                if self.q>7:
                    theta_f_guess = 150*np.pi/180.
            if self.chi_bh>=0.6 and self.chi_bh<0.7:
                if self.q<=5:
                    theta_f_guess = 30*np.pi/180. 
                if self.q>5 and self.q<=7:
                    theta_f_guess = 100*np.pi/180.
                if self.q>7:
                    theta_f_guess = 150*np.pi/180.
            if self.chi_bh>0.7 and self.chi_bh<0.9:
                if self.q<=5:
                    theta_f_guess = 86*np.pi/180
                if self.q>5:
                    theta_f_guess = 150*np.pi/180.
            if self.chi_bh>=0.9:
                if self.q<=4:
                    theta_f_guess = 85*np.pi/180
                if self.q>=4 and self.q<6:
                    theta_f_guess = 150*np.pi/180.
        if self.theta_i_deg == 180:
            chi_f = fsolve(self.ChifAligned, chi_f_guess, maxfev=10000)
            if chi_f >= 0:
                theta_f = np.pi
            if chi_f < 0:
                theta_f = 0
                chi_f = -chi_f
            m_f = self.M_f(chi_f,theta_f)
            return chi_f, theta_f*self.rad_deg, m_f, self.m_out
        
        chi_f, theta_f = fsolve(self.ChifThetaf, (chi_f_guess,theta_f_guess), maxfev=10000)
        m_f = float(self.M_f(chi_f,theta_f))
        
        return chi_f, theta_f*self.rad_deg, m_f, self.m_out

## OLD ##

    def secant(self,f):
        # SECANT METHOD
        # a and b are interval where I think the root of the function is included
        a = 0.01
        b = 0.99
        root = b - (b-a)*f(b)/(f(b)-f(a))
        it_sec = 1 # count the iterations of secant method
    
        while np.absolute(root-b)>self.x_err and np.absolute(f(root))>self.f_err and f(a)!=f(b) and it_sec<self.it_max:
            a = b
            b = root
            root = b - (b-a)/(f(b)-f(a))*f(b)
            it_sec = it_sec + 1
            #print(root)
        return root
