import matplotlib.pyplot as plt
import amplitude_module as am
import spin_module as sm
import ringdown.waveform as wf
import numpy as np

def PlotQNM(config):
    config_section_name='PLOT_QNM'
    qnm_par = {}
    qnm_par['l']=int(config[config_section_name]['l'])
    qnm_par['m']=int(config[config_section_name]['m'])
    qnm_par['n']=int(config[config_section_name]['n'])
    qnm_object = wf.QNM_fit(qnm_par['l'],qnm_par['m'],qnm_par['n'])
    qnm_par['x']=config[config_section_name]['x_var']
    qnm_par['y']=config[config_section_name]['y_var']
    qnm_par['m_bh_min']=float(config[config_section_name]['m_bh_min'])
    qnm_par['m_bh_max']=float(config[config_section_name]['m_bh_max'])
    m_bh=np.linspace(qnm_par['m_bh_min'], qnm_par['m_bh_max'], 9)
    qnm_par['chi_bh_min']=float(config[config_section_name]['chi_bh_min'])
    qnm_par['chi_bh_max']=float(config[config_section_name]['chi_bh_max'])
    chi_bh=np.linspace(qnm_par['chi_bh_min'], qnm_par['chi_bh_max'], 9)

    #plt.rcparams['xtick.top'] = plt.rcparams['xtick.labeltop'] = true 
    #plt.rcparams['ytick.right'] = plt.rcparams['ytick.labelright'] = true
    if qnm_par['x'] == "mass":
        for j in range(len(chi_bh)):
            y_var_arr = []
            y_var=0
            for k in range(len(m_bh)):
                if qnm_par['y'] == "freq":
                   y_var = qnm_object.f(m_bh[k],chi_bh[j])
                if qnm_par['y'] == "tau":
                   y_var = qnm_object.tau(m_bh[k],chi_bh[j])
                y_var_arr.append(y_var)
            plt.plot(m_bh, y_var_arr, label=r"$\chi_{BH}$: %.1f" %chi_bh[j])
        plt.xlabel(r"$M_{BH} [M_{\odot}]$")
        plt.ylabel(r"%s" %qnm_par['y'])
        plt.title("Mode = %d%d%d" %(qnm_par['l'], qnm_par['m'], qnm_par['n']), y=1.08)

    if qnm_par['x'] == "spin":
        for j in range(len(m_bh)):
            y_var_arr = []
            y_var=0
            for k in range(len(chi_bh)):
                if qnm_par['y'] == "freq":
                   y_var = qnm_object.f(m_bh[j],chi_bh[k])
                if qnm_par['y'] == "tau":
                   y_var = qnm_object.tau(m_bh[j],chi_bh[k])
                y_var_arr.append(y_var)
            plt.plot(m_bh, y_var_arr, label=r"$M_{BH}$: %d M_{odot} " %chi_bh[j])
        plt.xlabel(r"$\chi_{BH}$")
        plt.ylabel(r"%s" %qnm_par['y'])
        plt.title("Mode = %d%d%d" %(qnm_par['l'], qnm_par['m'], qnm_par['n']), y=1.08)
    plt.show()

def PlotAmplitude(config): 
    positions = [331,332,333,334,335,336,337,338]
    amplitudes = ["A22220", "A21210", "A33330", "A32320", "A32220", "A44440", "A43330", "A43430"]
    real = 0 # amplitude absolute value 
    #real = 1 # amplitude real part
    config_section_name='PLOT_AMP_QNM'
    qnm_par = {}
    qnm_par['m_ns']=float(config[config_section_name]['m_ns'])
    m_bh_min=float(config[config_section_name]['m_bh_min'])
    m_bh_max=float(config[config_section_name]['m_bh_max'])
    m_bh=np.linspace(m_bh_min, m_bh_max, 20)
    chi_bh_min=float(config[config_section_name]['chi_bh_min'])
    chi_bh_max=float(config[config_section_name]['chi_bh_max'])
    chi_bh=np.linspace(chi_bh_min, chi_bh_max, 9)

    for i in range(len(positions)):
        #plt.subplots_adjust(left=0.07, bottom=0.08, right=0.98, top=0.96, wspace=0.24, hspace=0.31)
        ax = plt.subplot(positions[i])
        for j in range(len(chi_bh)):
            qnm_par['chi_bh'] = chi_bh[j]
            amplitude = []
            for k in range(len(m_bh)):
                qnm_par['m_bh'] = m_bh[k]
                amplitude_object = am.Amplitude_class(qnm_par)
                switch={
                    "A22220": amplitude_object.A22220,
                    "A21210": amplitude_object.A21210,
                    "A33330": amplitude_object.A33330,
                    "A32320": amplitude_object.A32320, 
                    "A32220": amplitude_object.A32220,
                    "A44440": amplitude_object.A44440,
                    "A43330": amplitude_object.A43330,
                    "A43430": amplitude_object.A43430
                }
                a = switch[amplitudes[i]](real)
                amplitude.append(a)

            plt.plot(m_bh, amplitude, label=r"$\chi_{BH}$: %.2f" %chi_bh[j])
        plt.xlabel(r"$M_{BH}$")
        plt.ylabel("%s"%amplitudes[i])
        #chartBox = ax.get_position()
        #ax.set_position([chartBox.x0, chartBox.y0, chartBox.width, chartBox.height])
    plt.legend(loc='upper right', bbox_to_anchor=(1.6, 1), shadow=True, ncol=1)
    plt.show()
