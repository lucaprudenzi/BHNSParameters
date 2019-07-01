import matplotlib.pyplot as plt
import amplitude_module as am
import spin_module as sm
import numpy as np
import phenommixedamp_nospin as phenom_nospin
import phenommixedamp as phenom
import detectors as det
import lal
import ringdown.waveform as wf
from scipy.optimize import fsolve
def plotChiF_VarChi(config):
    config_section_name = 'CHI_F_VARCHI'
    binary_par = {} 
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 5*(q_max-q_min+1))
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 10*(chi_bh_max-chi_bh_min)+1)
    binary_par["type"]=config[config_section_name]['type']

    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    for j in range(len(chi_bh)):
        chi_f_arr = []
        binary_par['chi_bh'] = chi_bh[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            chi_f, theta_f, m_f, m_out = spin_object.FinalPar()
            chi_f_arr.append(chi_f)
        plt.plot(q, chi_f_arr, label=r"$\chi_{BH}:$ %.2f" %chi_bh[j])
        plt.legend()
    plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$=%d$^{\circ}$"%(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"]), y=1.08)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$\chi_f$")
    plt.show()


def plotChiF_VarTheta(config):
    config_section_name = 'CHI_F_VARTHETA'
    binary_par = {} 
    binary_par["chi_bh"]=float(config[config_section_name]['chi_bh'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 5*(q_max-q_min+1))
    theta_i_deg_min=float(config[config_section_name]['theta_i_deg_min'])
    theta_i_deg_max=float(config[config_section_name]['theta_i_deg_max'])
    theta_i_deg = np.linspace(theta_i_deg_min, theta_i_deg_max, 9)
    binary_par["type"]=config[config_section_name]['type']

    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    for j in range(len(theta_i_deg)):
        chi_f_arr = []
        binary_par['theta_i_deg'] = theta_i_deg[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            chi_f, theta_f, m_final, m_out = spin_object.FinalPar()
            chi_f_arr.append(chi_f)
        
        plt.plot(q, chi_f_arr, label= r'$\theta_i:$ %d$^{\circ}$' %theta_i_deg[j])
        plt.legend()
    plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\chi_{BH}$= %.2f "%(binary_par["m_ns"], binary_par["EOS"], binary_par["chi_bh"]), y=1.08)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$\chi_f$")
    plt.show()


def plotMOut(config):
    config_section_name = 'M_OUT'
    binary_par = {} 
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 5*(q_max-q_min+1))
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 10*(chi_bh_max-chi_bh_min)+1)
    binary_par["type"]=config[config_section_name]['type']

    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    for j in range(len(chi_bh)):
        m_out_arr = []
        binary_par['chi_bh'] = chi_bh[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            m_out = spin_object.M_out()
            m_out_arr.append(m_out)
        plt.plot(q, m_out_arr, label=r"$\chi_{BH}:$ %.2f" %chi_bh[j])
        plt.legend()
    plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$ = %d$^{\circ}$"%(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"]), y=1.08)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$M_{b,out}[M_{\odot}]$")
    plt.show()

def plotMOut2012(config):
    config_section_name = 'M_OUT_2012'
    binary_par = {} 
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 5*(q_max-q_min+1))
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 10*(chi_bh_max-chi_bh_min)+1)
    binary_par["type"]=config[config_section_name]['type']
    
    for j in range(len(chi_bh)):
        m_out_arr = []
        binary_par['chi_bh'] = chi_bh[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            m_out = spin_object.M_out2012()
            m_out_arr.append(m_out)
        plt.plot(q, m_out_arr, label=r"$\chi_{BH}:$ %.2f" %chi_bh[j])
        plt.legend()
    plt.title(r"2012 $M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$ = %d$^{\circ}$"%(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"]), y=1.08)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$M_{b,out}[M_{\odot}]$")
    plt.show()


def plotMf(config):
    config_section_name = 'M_F'
    binary_par = {} 
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 5*(q_max-q_min+1))
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 10*(chi_bh_max-chi_bh_min)+1)
    binary_par["type"]=config[config_section_name]['type']
   
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    for j in range(len(chi_bh)):
        m_f_arr = []
        binary_par['chi_bh'] = chi_bh[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            chi_f, theta_f, m_f, m_out = spin_object.FinalPar()
            m_f_arr.append(m_f/(binary_par["q"]*binary_par["m_ns"]+binary_par["m_ns"]))
        plt.plot(q, m_f_arr, label=r"$\chi_{BH}:$ %.2f" %chi_bh[j])
        plt.legend()
    plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$ = %d$^{\circ}$" %(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"]), y=1.08)

    plt.xlabel(r"$q$")
    plt.ylabel(r"$M_f$")
    plt.show()

def plotThetaF(config):
    config_section_name = 'THETA_F'
    binary_par = {} 
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 5*(q_max-q_min+1))
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 10*(chi_bh_max-chi_bh_min)+1)
    binary_par["type"]=config[config_section_name]['type']
    
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    for j in range(len(chi_bh)):
        theta_f_arr = []
        binary_par['chi_bh'] = chi_bh[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            chi_f, theta_f, m_f, m_out = spin_object.FinalPar()
            theta_f_arr.append(theta_f)
        plt.plot(q, np.abs(theta_f_arr), label=r"$\chi_{BH}:$ %.2f" %chi_bh[j])
        plt.legend()
    plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$ = %d$^{\circ}$"%(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"]), y=1.08)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$\theta_f$")
    plt.show()


def plotPhenomMix(config):
    config_section_name = 'PHENOM_MIX'
    binary_par = {} 
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["q"]=float(config[config_section_name]['q'])
    binary_par["chi_bh"]=float(config[config_section_name]['chi_bh'])
    binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    x_var = binary_par["x_var"] = config[config_section_name]['x_var']
    y_var = binary_par["y_var"] = config[config_section_name]['y_var']
    detector = config[config_section_name]['detector']
    dMpc = float(config[config_section_name]['distance'])
    
    spin_model = "spin"
    #spin_model = "nospin"
    binary_types = ["BBH", "BHNS"]

    G=6.67*10**(-8)
    c=3*10**10
    m_sol=2*10**33 
    parsec=3*10**18
    D = dMpc*10**6*parsec
    
    m_bh = binary_par["m_ns"]*binary_par["q"]
    m_ns = binary_par["m_ns"]
    M = m_bh+m_ns
    
    f_min = 10 #Hz
    f_max = 8000 #Hz
    f_arr = np.linspace(f_min, f_max, 1000)
    
    # frequency prefactor
    prefactor_f = G/c**3*m_sol*(m_ns+m_bh)
   
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True      
    ifo = det.Detector_class()
    if y_var == "fh" and x_var == "f":
        ifo.Plot(f_arr, detector)
    
    for binary_type in binary_types:
        amp_arr=[]
        binary_par["type"] = binary_type
        if binary_type == "BBH":
            color_type = "c"
        if binary_type == "BHNS":
            color_type = "r"
        
        if spin_model == "spin":
            phenom_object = phenom.PhenomBHNS_class(binary_par)
        if spin_model == "nospin":
            phenom_object = phenom_nospin.PhenomBHNS_nospin_class(binary_par)
        for j in range(len(f_arr)):
            amp = phenom_object.Amp_Phenom(f_arr[j])
            amp_abs=0
            if y_var == "fhD/m":   
                amp_abs = np.absolute(amp)*f_arr[j]*prefactor_f
            if y_var == "fh":
                amp_abs = np.absolute(amp)*f_arr[j]*prefactor_f*M*m_sol*G/c**2/D
            amp_arr.append(amp_abs)
        if x_var == "fm":
            plt.plot(prefactor_f*f_arr, amp_arr, label=binary_type, color=color_type)
        if x_var == "f":
            plt.plot(f_arr, amp_arr, label=binary_type, color=color_type)
        if y_var == "fh" and binary_type=="BHNS":
            snr = ifo.SNR(f_arr,amp_arr,detector) # print SNR
            print("SNR: %.2f" %snr)
    plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$ = %d$^{\circ}$, $q=$%d, $\chi_{BH}=$%.1f" %(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"], binary_par["q"], binary_par["chi_bh"]), y=1.08)
    
    if x_var == "fm":
        x_min = 0.005
        x_max = 0.10
        plt.xlim(x_min,x_max)
        plt.xlabel(r"$fm_0$")
    if y_var == "fhD/m":
        y_min = 0.01
        y_max = 0.2
        plt.ylim(y_min,y_max)
        plt.ylabel(r"$f|\tilde{h}(f)|D/m_0$")
    if x_var == "f":
        x_min = 50
        x_max = 8000
        plt.xlim(x_min, x_max)
        plt.xlabel(r"f (Hz)")
    if y_var == "fh":
        y_min = 0.001*M*m_sol*G/c**2/D
        y_max = 0.3*M*m_sol*G/c**2/D
        plt.ylim(y_min,y_max)
        plt.ylabel(r"$f|\tilde{h}(f)|$ at %d Mpc" %dMpc)
        
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()    

    plt.show()


def ParameterZones(config):
    config_section_name = 'PARAMETER_ZONES'
    binary_par = {}
    binary_par["theta_i_deg"]=float(config[config_section_name]['theta_i_deg'])
    binary_par["EOS"]=config[config_section_name]['EOS']
    binary_par["chi_ns"]=float(config[config_section_name]['chi_ns'])
    m_ns = binary_par["m_ns"]=float(config[config_section_name]['m_ns'])
    q_min=float(config[config_section_name]['q_min'])
    q_max=float(config[config_section_name]['q_max'])
    q = np.linspace(q_min, q_max, 20*(q_max-q_min+1))
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh = np.linspace(chi_bh_min, chi_bh_max, 15*(chi_bh_max-chi_bh_min)+1)
    binary_par["type"]=config[config_section_name]['type']
  
    G = 6.67*10**(-8)
    c = 3*10**(10)
    m_sol = 2*10**33

    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True 
    plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    m_out_chi = []
    m_out_q = []
    f_chi = []
    f_q = []
    qnm_object = wf.QNM_fit(2,2,0)

    for j in range(len(chi_bh)):
        m_out_arr = []
        m_out_first_time = 0
        f_tid_maj_f_220 = 0
        binary_par['chi_bh'] = chi_bh[j]
        for k in range(len(q)):
            binary_par['q'] = q[k] 
            spin_object = sm.Spin_class(binary_par)
            chi_f, theta_f, m_f, m_out = spin_object.FinalPar()
            m_out_arr.append(m_out)
            binary_par['q'] = q[k] 
            f_220 = qnm_object.f(m_f,chi_f)
            xi_tid = fsolve(spin_object.XiTidal,10)
            r_tid = xi_tid*(m_ns*m_sol*G/c**2)/spin_object.C*(1-2*spin_object.C)
            f_tid = float(1./(np.pi*(chi_f*m_f*m_sol*G/c**2+(r_tid**3/(m_f*m_sol*G/c**2))**(1./2)))*c)
            #print("")
            #print("m_out",m_out)
            #print("chi",chi_bh[j])
            #print("q",q[k])
            #print("")

            if f_tid>=f_220 and f_tid_maj_f_220==0:
                f_tid_maj_f_220=1
                f_chi.append(chi_bh[j])
                f_q.append(q[k])
                #print("")
                #print("f_220", f_220)
                #print("f_tid",f_tid)
                #print("chi",chi_bh[j])
                #print("q",q[k])
                #print("")
            if m_out==0 and m_out_first_time==0:
                m_out_first_time = 1
                m_out_chi.append(chi_bh[j])
                #print("")
                #print(chi_bh[j])
                m_out_q.append(q[k])
                #print(q[k])
                #print("")
        #plt.plot(q, m_out_arr, label=r"$\chi_{BH}:$ %.2f" %chi_bh[j])
        #plt.legend()
    #plt.title(r"$M_{NS}$ = %.2f $M_{\odot}$, EOS = %s, $\theta_i$ = %d$^{\circ}$" %(binary_par["m_ns"], binary_par["EOS"], binary_par["theta_i_deg"]), y=1.08)
    plt.plot(f_chi, f_q,label=r'$f_{tid}=f_{220}$')
    plt.plot(m_out_chi, m_out_q, label=r'$M_{OUT}=0$')
    plt.xlabel(r'$\chi_{BH}$')
    plt.ylabel(r'$q$')
    switch={
        "2H":[[0.1,9],[0.5,9],[0.4,6.5],[-0.01,6.0],[-0.01,5.8]],
        "H": [[0.2,8],[0.65,8],[0.6,4],[-0.05,3.4],[-0.05,3.2]],
        "HB":[[0.2,7],[0.6,6],[0.6,3.5],[-0.05,3],[-0.05,2.7]],
        "B": [[0.2,6],[0.65,5],[0.65,2.5],[-0.01,2.5],[-0.01,2.3]]
    }
    position = switch[binary_par["EOS"]]
    plt.text(position[0][0],position[0][1], r'$f_{220}<f_{tid}, M_{OUT}=0$')
    plt.text(position[1][0],position[1][1], r'$f_{220}<f_{tid}, M_{OUT}>0$')
    plt.text(position[2][0],position[2][1], r'$f_{tid}<f_{220}, M_{OUT}>0$')
    plt.text(position[3][0],position[3][1], r'$f_{tid}<f_{220}$')
    plt.text(position[4][0],position[4][1], r'$M_{OUT}=0$')
    plt.title(r"$M_{NS}=%.2f, EOS=%s, \theta_i=%d $"%(m_ns, binary_par["EOS"], binary_par["theta_i_deg"]),y=1.08)
    plt.legend()
    plt.show()
