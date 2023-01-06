import math                      # for mathematical operations
import numpy as npy              # for array procedures
import matplotlib.pyplot as plt  # plot
import os                        # to build path

IMF_List = ['Chabrier+03', 'Salpeter+55', 'Kroupa+93']

# Define IMF function
def IMF(m, imfRef):
    # IMF (per mass) : int_mlow^mup(IMF*dm) = 1.0
    # m must be in Msun
    if imfRef == 'Salpeter+55':
        A_Salpeter = 0.17163635690566686   # Intergartion between 0.1 and 100.
        return A_Salpeter*pow(m,-1.35)
    elif imfRef == 'Chabrier+03':
        A_Chabrier = 0.8465815581229609  #0.85   # for m_low = 0.1
        B_Chabrier = 0.23903479288177717 #0.24
        m_c        = 0.079  # in Msun
        s          = 0.69
        if m <= 1.0:    # in Msun
            return A_Chabrier*math.e**(-(math.log10(m)-math.log10(m_c))**2./(2.*s**2.))
        elif m > 1.0:  # in Msun
            return B_Chabrier*m**(-1.3)
    elif imfRef == 'Kroupa+93':
        A_Kroupa = 0.5797490548994692 #0.58
        B_Kroupa = 0.30986587417040595 #0.31
        C_Kroupa = 0.30986587417040595 #0.31
        if m < 0.5:  # in Msun
            return A_Kroupa*pow(m,-0.3)
        elif m < 1.0: # in Msun
            return B_Kroupa*pow(m,-1.2)
        elif m >= 1.0: # in Msun
            return C_Kroupa*pow(m,-1.7)
    elif imfRef == 'Scalo+98':
        A_Scalo = 0.3941196163613269  #0.39
        B_Scalo = 0.3941196163613269  #0.39
        C_Scalo = 0.16169009902003156 #0.16
        if m < 1.0: # in Msun
            return A_Scalo*m**(-0.2)
        elif m < 10.0: # in Msun
            return B_Scalo*m**(-1.7)
        elif m >= 10.0: # in Msun
            return C_Scalo*m**(-1.3)

# Test IMF
# create a Mass scale and test int_Mmin^Mmax(Phi(m))
def test():
    # define IMF list
    IMF_list = ['Salpeter+55', 'Chabrier+03', 'Kroupa+93', 'Scalo+98']
    #
    # define mass scale
    # We consider the following mass range:  M in [0.1 : 100] Msun
    M_min = 0.1
    M_max = 100.0
    dm = 0.0001
    Nsteps = int(math.floor((M_max - M_min)/dm) + 1)
    Mscale = npy.linspace(M_min, M_max, Nsteps)
    #
    for l, item in enumerate(IMF_list):
        print('--> imfRef: %s' %(IMF_list[l]))
        i = 0.
        for m in range(len(Mscale) - 1):
            Md = Mscale[m]
            int_d = IMF(Md, imfRef=IMF_list[l])
            #
            Mu = Mscale[m + 1]
            int_u = IMF(Mu, imfRef=IMF_list[l])
            #
            dM = abs(Mu - Md)
            i += 0.5 * dM * (int_d + int_u)
        print(' --> int = %.15f' %(i))
        
# plot 
def plot():
    #
    # Create a mass [log]-scale from 0.1 to 100.
    M_min = 0.1
    M_max = 100.0
    Mscale = 10.**(npy.linspace(math.log10(M_min), math.log10(M_max), 1000))
    IMF_W = []
    IMF_list = ['Salpeter+55','Chabrier+03','Kroupa+93','Scalo+98']
    #
    # plot
    plt.figure(figsize = (5.5, 5))
    color = ['blue', 'green', 'orange', 'red']  # color reference for different IMF
    for i, item in enumerate(IMF_list):
        IMF_W.append([])  # Add a new IMF
        for m in Mscale:
            IMF_W[i].append(IMF(Mscale[m], imfRef = IMF_list[i]))
        plt.plot(Mscale,IMF_W[i], c=color[i], label=IMF_list[i])
    plt.xscale('log')           # use log scale for x axis
    plt.yscale('log')           # use log scale for y axis
    plt.xlabel('$m_{\star}$ ($M_{\odot}$)')   # label for x axis
    plt.ylabel('$\phi(m_{\star})$')           # label for y axis
    plt.grid()
    plt.legend()
    # Save the plot in pdf format
    # path eGalICS
    plot_path = 'sp/plots'
    filename = 'IMFs.pdf'
    filename = os.path.join(plot_path, filename)
    plt.savefig(filename)

if __name__ == "__main__":
    test()
    plot()