import numpy as np
import ringdown.waveform as wf
import spin_module as sm
import matplotlib.pyplot as plt
import lal

#From arxiv:gr-qc/0512160v2 eq. 4.12a
def SigmaChiF_Fisher(qnm, f_rd, tau_rd, chi, snr):
    Q = np.pi*f_rd*tau_rd
    sigma = 1./snr*abs(2*Q/(qnm.q_coeff[1]*qnm.q_coeff[2]*\
                         (1-chi)**(qnm.q_coeff[2]-1))*(1+1./(16*Q**2)))
    return sigma

#From arxiv:gr-qc/0512160v2 eq. 4.12b
def SigmaMF_Fisher(qnm, f_rd, tau_rd, chi_f, m_f, snr):
    Q = np.pi*f_rd*tau_rd
    fprime = lal.C_SI**3/(2.*np.pi*lal.G_SI*m_f*lal.MSUN_SI)*\
            (-qnm.f_coeff[1]*qnm.f_coeff[2]*(1-chi_f)**(qnm.f_coeff[2]-1))
    Qprime = -qnm.q_coeff[1]*qnm.q_coeff[2]*(1-chi_f)**(qnm.q_coeff[2]-1)
    sigma = 1./snr*abs(2*m_f*Q*fprime/(f_rd*Qprime)*(1+1./(16*Q**2)))
    return sigma

# Keeping chi_f fixed, compare the error on the final mass
# with the mass of the disk produced by the same system

def sigmaMF(config):
    colors = ['b','g','r','c','m','y','k','hotpink','olivedrab','gray']
    sigma_par={}
    config_section_name='SIGMA_MF' 
    sigma_par['m_ns']=float(config[config_section_name]['m_ns'])
    sigma_par['type']=config[config_section_name]['type']
    sigma_par['theta_i_deg']=float(config[config_section_name]['theta_i_deg'])
    sigma_par['EOS']=config[config_section_name]['EOS']
    snr=float(config[config_section_name]['snr'])
    l=float(config[config_section_name]['l'])
    m=float(config[config_section_name]['m'])
    n=float(config[config_section_name]['n'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q_arr = np.linspace(q_min,q_max,30)
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 10)
    qnm = wf.QNM_fit(l,m,n)
    
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True      
    for i in range(len(chi_bh)):
        sigma_par['chi_bh']=chi_bh[i]
        m_out_arr = []
        err_m_f = []
        for j in range(len(q_arr)):
            sigma_par['q'] = q_arr[j]
            spin_object = sm.Spin_class(sigma_par)
            chi_f, theta_f, m_f, m_out = spin_object.FinalPar() 
            f_rd = qnm.f(m_f, chi_f)
            tau_rd = qnm.tau(m_f, chi_f)
            sigma_m_f = SigmaMF_Fisher(qnm, f_rd, tau_rd, chi_f, m_f, snr)
            m_out_arr.append(m_out)
            err_m_f.append(sigma_m_f)
        plt.plot(q_arr, err_m_f, color=colors[i%10], label=r"$\sigma_M$, $\chi_{BH}$:%.1f"%chi_bh[i], linestyle='--')
        plt.plot(q_arr, m_out_arr, colors[i%10], label=r"$M_{OUT}$,$\chi_{BH}:%.1f$"%chi_bh[i], linestyle='-')
        #idx = np.argwhere(np.diff(np.sign(np.array(err_m_f) - np.array(m_out_arr))).flatten()
        #plt.plot(q_arr[idx], m_out_arr[idx], 'ro')
    plt.title(r"$M_{NS}$=%.2f, EOS= %s, mode = %d%d%d, SNR = %d"%(sigma_par["m_ns"], sigma_par["EOS"], l,m,n, snr), y=1.08)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$[M_{\odot}]$")
    plt.legend()
    plt.show()


########################
## DEPRECATED METHODS ##
########################
def sigma_freq(snr, mode1, mode2, binary_par):
    
    spin_object = sm.Spin_class(binary_par)
    a_f, theta_f, m_f, m_out = spin_object.FinalPar()
    
    if mode1 == mode2:
        print("non orthogonal mode")
    else:
        l1 = int(mode1[0]); m1 = int(mode1[1]); n1 = int(mode1[2])
        l2 = int(mode2[0]); m2 = int(mode2[1]); n2 = int(mode2[2])
    
    # Pair of modes
    qnm1 = wf.QNM_fit(l1,m1,n1)
    qnm2 = wf.QNM_fit(l2,m2,n2)
    
    f1 = qnm1.f(m_f, a_f)
    f2 = qnm2.f(m_f, a_f)
    
    tau1 = qnm1.tau(m_f, a_f)
    tau2 = qnm2.tau(m_f, a_f)

    Q1 = np.pi*f1*tau1
    Q2 = np.pi*f2*tau2

    amplitude_object = am.Amplitude_class(binary_par)
    switch={
        "220": amplitude_object.A22220,
        "210": amplitude_object.A21210,
    }
    
    A1 = switch[mode1](0)
    A2 = switch[mode2](0)
    #print(f1,f2,tau1,tau2,A1,A2, m_f, a_f)
    sigma = 1./(snr*2*2**(1./2))*(f1**3*(3+16*Q1**4)/(A1**2*Q1**7)*\
                                   (A1**2*Q1**3/(f1*(1+4*Q1**2))+\
                                    A2**2*Q2**3/(f2*(1+4*Q2**2))))**(1./2)
    sigma = 1./(2**(1./2)*np.pi*tau1*snr)
    return sigma

def sigma_tau(snr, mode1, mode2, binary_par):
    
    spin_object = sm.Spin_class(binary_par)
    a_f, theta_f, m_f, m_out = spin_object.FinalPar()
    
    if mode1==mode2:
        print("Error, non orthogonal modes")
    else:
        l1 = int(mode1[0]); m1 = int(mode1[1]); n1 = int(mode1[2])
        l2 = int(mode2[0]); m2 = int(mode2[1]); n2 = int(mode2[2])
    
    qnm1 = wf.QNM_fit(l1,m1,n1)
    qnm2 = wf.QNM_fit(l2,m2,n2)
    
    f1 = qnm1.f(m_f, a_f)
    f2 = qnm2.f(m_f, a_f)
    
    tau1 = qnm1.tau(m_f, a_f)
    tau2 = qnm2.tau(m_f, a_f)

    Q1 = np.pi*f1*tau1
    Q2 = np.pi*f2*tau2

    amplitude_object = am.Amplitude_class(binary_par)
    switch={
        "220": amplitude_object.A22220,
        "210": amplitude_object.A21210,
    }
    
    A1 = switch[mode1](0)
    A2 = switch[mode2](0)
    sigma = 2./(snr*np.pi)*((3+4*Q1**2)/(A1**2*f1*Q1)*\
                            (A1**2*Q1**3/(f1*(1+4*Q1**2))+\
                             A2**2*Q2**3/(f2*(1+4*Q2**2))))**(1./2)
    sigma = 2.*tau1/snr
    
    return sigma

def dQdTau(f):
    dqdtau = np.pi*f
    return dqdtau

def dQdF(tau):
    dqdf = np.pi*tau
    return dqdf

def sigma_Q(f, sigma_f, tau, sigma_tau):
    dqdtau = dQdTau(f)
    dqdf = dQdF(tau)
    sigma_q = (dqdf**2*sigma_f**2+dqdtau**2*sigma_tau**2)
    return sigma_q

def dAdQ(qnm, Q):
    dadq = -1./(qnm.q_coeff[1]*qnm.q_coeff[2])*\
            ((Q-qnm.q_coeff[0])/qnm.q_coeff[1])**(1./qnm.q_coeff[2]-1)
    return dadq

def sigma_A(qnm, Q, sigma_Q):
    dadq = dAdQ(qnm, Q)
    sigma_a = (dadq**2*sigma_Q**2)**(1./2)
    return sigma_a

def dMdF(qnm, f, a):
    prefactor = lal.C_SI*lal.C_SI*lal.C_SI/(2.*np.pi*lal.G_SI*lal.MSUN_SI)
    dmdf = prefactor*(-1./f**2*(qnm.f_coeff[0]+qnm.f_coeff[1]*\
                                (1-a)**qnm.f_coeff[2]))
    return dmdf

def dMdA(qnm, f, a):
    prefactor = lal.C_SI*lal.C_SI*lal.C_SI/(2.*np.pi*lal.G_SI*lal.MSUN_SI)
    dmda = prefactor*(-1./f*(qnm.f_coeff[1]*qnm.f_coeff[2]*\
                             (1-a)**(qnm.f_coeff[2]-1)))
    return dmda

#sigma_m_final = (dmdf**2*sigma_f**2+dmda**2*sigma_a**2+2*sigma_a*sigma_f*dmdf*dmda)**(1./2)
def Sigma_M_Final(qnm, f, sigma_f, tau, sigma_tau, a):
    dmdf = dMdF(qnm, f, a)
    Q = np.pi*f*tau
    sigma_q = sigma_Q(f, sigma_f, tau, sigma_tau)
    dmda = dMdA(qnm, f, a)
    sigma_a = sigma_A(qnm, Q, sigma_q)
    sigma_m_final = (dmdf**2*sigma_f**2+dmda**2*sigma_a**2)**(1./2)

    return sigma_m_final


