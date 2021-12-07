#############################################################
# For given starting parameters, this code finds the optimal
# and mu_CV for a specified Eg and plots the JV curve and the
# efficiency vs mu_CV.
#
# Author: Jordan Pagé (EN/FR; he/him)
# Email:  jpage019@uottawa.ca
# Date:   December 7th 2021
#############################################################

import numpy as np
import scipy
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from scipy import optimize

def main():
    # Parameters
    V = 0               # Starting voltage across cell
    q = 1.60217e-19     # electron charge
    mu_CV = V           # Chemical potential of photon

    Tc = 300            # K, tempurature of cell
    Ts = 6000           # K, tempurature of sun surface

    eG = 1              # energy of gap, 1 eV
    eM = 20             # 20 eV
    em = eG             # lower bound of integration

    O_sb = 5.670374419e-8 # Stefan-Boltzmann constant, W/(m²K⁴)

    # creating array of mu for plotting the general look of the function
    lst = []
    step = 0.2 # The function has a very steep slope near e = mu_CV, so I want more points
               # taken the closer we get to that point
    while mu_CV<em*0.99982:
        lst.append(mu_CV)
        mu_CV+=(eG-mu_CV)*step  # updating mu


    J_v = Jv(lst, em, eM, Ts = Ts, Tc = Tc)

    maximum = [optimize.fminbound(neg_Jv, 0, eG, args = (em, eM))]
    #maximum = optimize.fmin(neg_Jv, [0.9], args = (em, eM))

    maximum = optimize.minimize(neg_Jv, [0.9], args = (em, eM), bounds = [(0,1)], tol = 1e-9, options={'ftol':1e-9,'gtol':1e-9} )
    print('')
    print(maximum)
    """
  message: 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
     nfev: 6
      nit: 1
     njev: 3
   status: 0
  success: True
    """
    maximum = maximum.x

    ideal_mu = Jv(maximum, em, eM)
    print(maximum)
    print(ideal_mu)

    efficiency = [i/(O_sb*Ts**4) for i in J_v] # Pretty sure this doesn' give the correct result...


    # Plotting power curve
    fig, axs = plt.subplots(1)
    plt.axes(axs)
    plt.plot(lst, J_v, c='blue', linewidth = 2, label = r'Output Power')
    plt.scatter(maximum, ideal_mu, c='b', label = f"({maximum[0]:.5}; {ideal_mu[0]:.5})")
    plt.xlabel(r'$\mu_{CV}~(eV)$')
    #plt.ylabel(r'Photon Flux $eV/(eV C m²s)$')
    plt.legend(fontsize = 12)
    plt.grid()
    plt.show()

    # Plotting efficiency
    fig, axs = plt.subplots(1)
    plt.axes(axs)
    plt.plot(lst, efficiency, c='blue', linewidth = 2, label = r'Efficiency')
    #plt.scatter(maximum, ideal_mu, c='b', label = f"({maximum:.5}; {ideal_mu:.5})")
    #plt.xlabel(r'$\mu_{CV}~(eV)$')
    #plt.ylabel(r'Photon Flux $eV/(eV C m²s)$')
    plt.legend(fontsize = 12)
    plt.grid()
    plt.show()


def N_dot(e, mu_CV, T):
    # Some useful constants
    pi = np.pi
    h = 4.135667696e-15 # Planck's constant in eV s
    c = 299792458       # Speed of ligh in a vacuum, m/s
    kB = 8.617333262145e-5 # Boltzmann's constant in eV/K
    q = 1.60217e-19

    # Defining efficiency function
    C = 2*pi/(h**3 * c**2) # 1/((eV s)³ (m/s)²) = s²/(eV³s³m²) = 1/eV³m²s
    return C*(e)**2 / (np.exp((e-mu_CV)/(kB*T))-1) # Photon flux for sun
    # in the exponent: eV/(K*eV/K) = 0
    # (1/eV³m²s) * (eV)² = 1/(eV m²s)


def Iq(mu_CV, em, eM, Ts = 6000, Tc = 300):
    q = 1.60217e-19
    integral1 = integrate.quad(N_dot, em, eM, args = (0, Ts))[0] # Flux-in from Sun
    integral2 = integrate.quad(N_dot, em, eM, args = (mu_CV, Tc))[0] # Flux-out
    return (integral1 - integral2)*q

def Jv(mu_CV_lst, em, eM, Ts = 6000, Tc = 300):
    if isinstance(mu_CV_lst, list) or isinstance(mu_CV_lst, np.ndarray):
        J_v = [Iq(i, em, eM)*i for i in mu_CV_lst] # This line is used for plotting a curve
    else:
        J_v = Iq(mu_CV_lst, em, eM)*mu_CV_lst # This line is used when optimize.fminbound() calls the function
    return J_v

def neg_Jv(mu_CV_lst, em, eM, Ts = 6000, Tc = 300):
    # This function is the negative of Jv(). It is used by fminboun() to find the max of Jv()
    if isinstance(mu_CV_lst, list) or isinstance(mu_CV_lst, np.ndarray):
        #J_v = [-Iq(i, em, eM)*i for i in mu_CV_lst] # This line is used for plotting a curve
        J_v = -Iq(mu_CV_lst, em, eM)*mu_CV_lst
    else:
        J_v = -Iq(mu_CV_lst, em, eM)*mu_CV_lst # This line is used when optimize.fminbound() calls the function
    return J_v

main()