############################################
# Detailed balance calculations for IBSC
#
# Author: Jordan Pagé
# email: jpage019@uottawa.ca
# January 7th 2022
############################################

import numpy as np
import scipy
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from scipy import optimize
import csv

#Ideal mu_CV: 1.8661344489992717
#Ideal eG:    1.9886288598286463
#Ideal eI:    0.7286495007361697
#Ideal eC:    1.2599793590924766
#Efficieny:   63.16%
#Itterations: 1916

def main():
    # Parameters
    sigma_sb = 5.670374419e-8 # Stefan-Boltzmann constant, W/(m²K⁴)
    Tc = 300   # K, tempurature of cell
    Ts = 6000  # K, tempurature of sun surface

    # all units are in eV
    initial_eG = 2.                 # initial guess of optimal energy of gap
    initial_eI = initial_eG/4       # initial guess of optimal energy of intermediate band
    eM = 20                         # "infinity", 20eV

    #mu_CV = 1.8661344489992717
    #eG =    1.9886288598286463
    #eI =    0.7286495007361697
    #eC =    1.2599793590924766
    print('Options:\n 1. Find global maxima\n 2. Reproduce Fig. 2 from the Luque & Martì paper\n 3. Perform both options 1 & 2\n 4. Exit\n')
    option = str(input('Enter "1", "2" or "3" \n>>> '))
    
    if option == '4':
        return

    if option != '1' and option != '2' and option != '3':
        print(f'{option} is not a valid input. Program terminating now.')

    if option == '1' or option == '3':
        # Finding optimal eG, eI & mu_CV
        ideal_eG, ideal_mu_CV, ideal_eI = find_globalPeak(initial_eG, initial_eI, eM)
        ideal_eC = ideal_eG - ideal_eI

        # creating a grid for plotting
        lst = np.linspace(0,ideal_eG*0.99,200)

        # making curve
        #params = [ideal_eG, lst, ideal_eI]
        power = P(ideal_eG, lst, ideal_eI, eM, Ts = Ts, Tc = Tc)

        # finding max power
        #ideal_mu_CV = optimise_muCV(ideal_eG, ideal_eI, eM) # mu_CV corresponding to max power

        #params = [ideal_eG, ideal_mu_CV, ideal_eI]
        maximum = P(ideal_eG, ideal_mu_CV, ideal_eI, eM, Ts = Ts, Tc = Tc) # max power

        eff = [i/(sigma_sb*Ts**4) for i in power] # max efficiency
        peak_eff = maximum/(sigma_sb*Ts**4) # max efficiency
        ######################## Plotting  #########################
        # Plotting power curve

        fig, axs = plt.subplots(1)
        plt.axes(axs)
        plt.plot(lst, power, c='blue', linewidth = 2, label = r'Output current')
        plt.scatter(ideal_mu_CV, maximum, c='b', label = f"({ideal_mu_CV:.5}; {maximum:.5})")
        plt.xlabel(r'$\mu_{CV}~(eV)$')
        plt.title(r'$\epsilon_g=$'+f'{ideal_eG:.3} eV,  ' r'$\epsilon_i=$' +f'{ideal_eI:.3} eV, '+ r'$\mu_{CV}=$'+f'{ideal_mu_CV:.3} eV')
        #plt.ylabel(r'Photon Flux $eV/(eV C m²s)$')
        plt.legend(fontsize = 12)
        plt.grid()
        plt.show()

        print('')
        print(f'Ideal mu_CV: {ideal_mu_CV}')
        print(f'Ideal eG:    {ideal_eG}')
        print(f'Ideal eI:    {ideal_eI}')
        print(f'Ideal eC:    {ideal_eC}')


        print('Efficieny:   ' + str(round(peak_eff*100,2)) + '%')


        # Plotting efficiency
        fig, axs = plt.subplots(1)
        plt.axes(axs)
        plt.plot(lst, eff, c='blue', linewidth = 2, label = r'Efficiency')
        plt.title(r'Efficiency of an IB solar cell with $\epsilon_C$ ='+f' {round(ideal_eC,3)}'+r'; $\epsilon_G = $'+f'{round(ideal_eG,3)}')
        #plt.scatter(maximum, ideal_mu, c='b', label = f"({maximum:.5}; {ideal_mu:.5})")
        plt.xlabel(r'$\mu_{CV}~(eV)$')
        plt.ylabel(r'Efficiency $\eta$')
        #plt.legend(fontsize = 12)
        plt.grid()
        plt.show()

    if option == '2' or option == '3':

        def prob_constraints2(eG, eI, eM):
            # constraints eG and eI must have
            return 0.49*eG - eI # this constrains eI to the range [0, 0.49*eG]. When eI gets too close to eG/2, things get weird.

        ### Finding optimal eG & eI for given eI (this is to recreate fig 2 from Luque & Marti)
        eIs = np.linspace(0.5, 1.3, 9) # the values of eI explicitly plotted in the Luque & Marti paper fig 2
        initial_eGs = eIs * 2.1        # the starting guess values of eG
        opt_eGs = [optimize.minimize(optimise_eG, i, args = (j, eM,), constraints = {'type': 'ineq', 'fun': prob_constraints2, 'args': (j, eM,)}, method = 'COBYLA', tol = 1e-9 ).x for i,j in zip(initial_eGs, eIs) ]
        opt_mu_CV = [optimise_muCV(i, j, eM) for i,j in zip(opt_eGs, eIs)]

        #print('        Ei:', eIs)
        #print('        Eg:', opt_eGs)
        opt_eGs = [float(i) for i in opt_eGs]
        #print('        Eg:', opt_eGs)
        #print('Initial Eg:', initial_eGs)

        power = P2(opt_eGs, opt_mu_CV, eIs, eM)
        eff = [i/(sigma_sb*Ts**4) for i in power] # max efficiency
        #print('     Power:', power)
        #print('Efficiency:', eff)

        """
        # Plotting power curve with changing eG
        fig, axs = plt.subplots(1)
        plt.axes(axs)
        plt.plot(eGs, efficiency, c='blue', linewidth = 2, label = r'Peak efficiency for given $\epsilon_G$')
        plt.scatter(ideal_eG, peak_eff, c='b', label = f"Global peak efficiency: ({ideal_eG:.5}; {peak_eff:.5})")
        plt.xlabel(r'$\epsilon_{G}~(eV)$')
        plt.ylabel(r'Efficiency $\eta$')
        plt.legend(fontsize = 12)
        plt.grid()
        plt.show()
        """

        #### LUQUE & MARTI fig 2 ####
        # Plotting power curve
        fig, axs = plt.subplots(1)
        plt.axes(axs)
        plt.plot(eIs, power, c='blue', linewidth = 2, label = r'Output power')
        #plt.scatter(ideal_mu_CV, maximum, c='b', label = f"({ideal_mu_CV:.5}; {maximum:.5})")
        plt.xlabel(r'$\epsilon_i$')
        #plt.title(r'$\epsilon_g=$'+f'{ideal_eG:.3} eV,  ' r'$\epsilon_i=$' +f'{ideal_eI:.3} eV, '+ r'$\mu_{CV}=$'+f'{ideal_mu_CV:.3} eV')
        #plt.ylabel(r'Photon Flux $eV/(eV C m²s)$')
        plt.legend(fontsize = 12)
        plt.grid()
        plt.show()

        # Plotting power curve
        fig, axs = plt.subplots(1)
        plt.axes(axs)
        plt.plot(eIs, eff, c='blue', linewidth = 2)#, label = r'Efficiency')
        for i in range(0,len(opt_eGs)):
            axs.annotate( str(round(opt_eGs[i],2)), xy = (eIs[i], eff[i]) )
        #plt.scatter(ideal_mu_CV, maximum, c='b', label = f"({ideal_mu_CV:.5}; {maximum:.5})")
        plt.xlabel(r'$\epsilon_i$')
        #plt.title(r'$\epsilon_g=$'+f'{ideal_eG:.3} eV,  ' r'$\epsilon_i=$' +f'{ideal_eI:.3} eV, '+ r'$\mu_{CV}=$'+f'{ideal_mu_CV:.3} eV')
        plt.ylabel('Efficiency')
        #plt.legend(fontsize = 12)
        plt.savefig(f'Luque&Marti_Fig2_efficiency.png')
        plt.grid()
        plt.show()

        save = input('Would you like to save the data? (y/n) >>> ')
        while save != 'y' and save != 'n':
            print('That is not valid input.')
            save = input('Would you like to save the data? (y/n) >>> ')

        if save == 'y':
            header = ['Eg', 'Ei', 'mu_CV', 'Power', 'Efficiency']
            with open(f'Luque&Marti_Fig2.txt', 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=header)
                writer.writeheader()
                for i,g,m,p,e in zip(eIs, opt_eGs, opt_mu_CV, power, eff ):
                    print(round(i,3),round(g,3),round(m,3),round(p,3),round(e,3))
                    writer.writerow({'Eg': str(g), 'Ei': str(i), 'mu_CV': str(m), 'Power': str(p), 'Efficiency': str(e)})

    ############################################################

    # when testing if mu_CV optimisation works, comment out the section bellow (to line 50) and uncomment lines 52-54
def find_globalPeak(initial_eG, initial_eI, eM, Ts = 6000, Tc = 300):
    ### Method 3: using using 1D & 2D optimizer ###
    def prob_constraints2(params, eM):
        # constraints eG and eI must have
        eG = params[0]
        eI = params[1]
        return 0.49*eG - eI # this constrains eI to the range [0, 0.49*eG]. When eI gets too close to eG/2, things get weird.

    cons = {'type': 'ineq',
             'fun': prob_constraints2,
             'args': (eM,)}

    # finding the optimal values of eG and eI
    params = [initial_eG, initial_eI]
    opt_params = optimize.minimize(optimise_all, params, args = (eM), constraints = cons, method = 'COBYLA', options = {'maxiter':2000,} )

    print(opt_params)
    print('')

    ideal_eG = opt_params.x[0]
    ideal_eI = opt_params.x[1]
    ideal_eC = ideal_eG - ideal_eI

    # finding the optimal mu_CV for the now ideal eG and eI values
    ideal_mu_CV = optimise_muCV(ideal_eG, ideal_eI, eM)
    #ideal_mu_CV = 2.145
    # when testing if mu_CV optimisation works, comment out the section above ( to line 23) and uncomment lines 52-54
    return ideal_eG, ideal_mu_CV, ideal_eI


############## Mapping mu_CI to arctan(x) ################
# the range of mu_CI is [mu_CV-eI, eG-eI=eC]
# the range of y = arctan(x) is [-inf, inf] -> [-pi/2, pi/2]
#              y = arctan(x)+pi/2 -> [-inf, inf] -> [0, pi]
#              has a "width" of pi
#
# "width" of mu_CI = max-min = (eG-eI)-(mu_CV-eI) = eG-mu_CV
# need to map pi -> mu_CV-eG:
# (arctan(x)+pi/2)/pi -> this goes form 0 to 1
# (eG-mu_CV)*(arctan(x)+pi/2)/pi -> this goes from 0 to (eG-mu_CV)
# We want it to go from mu_CV-eI to eG-eI, so be need to shift it by mu_CV-eI
# mu_CI = (eG-mu_CV)*(arctan(x)+pi/2)/pi + mu_CV-eI
# -> x = tan( (mu_CI-(mu_CV-eI)) * pi/(eG-mu_CV) - pi/2 )

def x_to_mu(x, mu_CV, eG, eI):
    x = (eG-mu_CV)*(np.arctan(x) + np.pi/2)/np.pi + mu_CV-eI
    return x

def mu_to_x(mu, mu_CV, eG, eI):
    mu = np.tan( (mu-(mu_CV-eI)) * (np.pi/(eG-mu_CV)) - np.pi/2 )
    return mu

##################### photon flux ########################

def N_dot(e, T, mu):
    # Some useful constants
    pi = np.pi
    h = 4.135667696e-15 # Planck's constant in eV s
    c = 299792458       # Speed of ligh in a vacuum, m/s
    kB = 8.617333262145e-5 # Boltzmann's constant in eV/K
    q = 1.60217e-19     # charge of electron

    # Defining efficiency function
    C = 2*pi/(h**3 * c**2) # 1/((eV s)³ (m/s)²) = s²/(eV³s³m²) = 1/eV³m²s
    return C*(e)**2 / (np.exp((e-mu)/(kB*T))-1) # Photon flux for sun
    # in the exponent: eV/(K*eV/K) = 0
    # (1/eV³m²s) * (eV)² = 1/(eV m²s)

################## Solving for mu_CI #####################

def eq7(x, eG, mu_CV, eC, eI, Ts = 6000, Tc = 300):
    # This is based off of eq 7 in the Luque & Marti paper

    # unmapping mu_CI to use it in the integrals
    mu_CI = x_to_mu(x, mu_CV, eG, eI) #original line
    #mu_CI = mu_to_x(x, mu_CV, eG, eI)

    mu_IV = mu_CV - mu_CI # eq 8 of Luque & Marti

    integral1 = integrate.quad(N_dot, eI, eC, args = (Ts, 0))[0] # Flux-in from Sun
    integral2 = integrate.quad(N_dot, eI, eC, args = (Tc, mu_IV))[0] # Flux-out from V->C
    integral3 = integrate.quad(N_dot, eC, eG, args = (Ts, 0))[0]
    integral4 = integrate.quad(N_dot, eC, eG, args = (Tc, mu_CI))[0] # Flux-out from IB->C

    return integral1 - integral2 - (integral3 - integral4)

################ eq 6 from Luque & Marti #################

def Iq(mu_CV, eG, eC, eI, eM, Ts = 6000, Tc = 300):
    q = 1.60217e-19

    # mapping mu_CI to a arctan() function
    mu_CI_initial_guess = (mu_CV-eI)/2
    x_initial_guess = mu_to_x(mu_CI_initial_guess, mu_CV, eG, eI) #original line
    #x_initial_guess = x_to_mu(mu_CI_initial_guess, mu_CV, eG, eI)

    # solving for the mapped vaeibale
    ideal_x = optimize.fsolve( eq7, x_initial_guess, args = ( eG, mu_CV, eC, eI ) )

    # unmapping mu_CI to a arctan() function
    mu_CI = x_to_mu(ideal_x, mu_CV, eG, eI) # original line
    #mu_CI = mu_to_x(ideal_x, mu_CV, eG, eI)

    mu_IV = mu_CV - mu_CI # eq 8 of Luque & Marti

    integral1 = integrate.quad(N_dot, eG, eM, args = (Ts, 0))[0] # Flux-in from Sun
    integral2 = integrate.quad(N_dot, eG, eM, args = (Tc, mu_CV))[0] # Flux-out from V->C
    integral3 = integrate.quad(N_dot, eC, eG, args = (Ts, 0))[0]
    integral4 = integrate.quad(N_dot, eC, eG, args = (Tc, mu_CI))[0] # Flux-out from IB->C

    return (integral1 - integral2 + (integral3 - integral4))*q

####################### Power curve ######################

def P(eG, mu_CV, eI, eM, Ts = 6000, Tc = 300):
    eC = eG - eI

    if isinstance(mu_CV, list) or isinstance(mu_CV, np.ndarray):
        power = [Iq(i, eG, eC, eI, eM)*i for i in mu_CV] # This line is used for plotting a curve
    else:
        power = Iq(mu_CV, eG, eC, eI, eM)*mu_CV # This line is used when optimize.fminbound() calls the function
    return power

def P2(eG, mu_CV, eI, eM, Ts = 6000, Tc = 300):
    eC = eG - eI

    if isinstance(mu_CV, list) or isinstance(mu_CV, np.ndarray):
        power = [Iq(i, j, k, l, eM)*i for i,j,k,l in zip(mu_CV, eG, eC, eI)] # This line is used for plotting a curve
    else:
        power = Iq(mu_CV, eG, eC, eI, eM)*mu_CV # This line is used when optimize.fminbound() calls the function
    return power

################### Negative power curve ###################

def neg_P_mu(mu_CV, eG, eI, eM, Ts = 6000, Tc = 300):
    # This function is the negative of Jv(). It is used by fminboun() to find the max of Jv()
    # I tried "return -P()" but Python it didn't like that

    eC = eG - eI

    if isinstance(mu_CV, list) or isinstance(mu_CV, np.ndarray):
        power = [-Iq(i, eG, eC, eI, eM)*i for i in mu_CV] # This line is used for plotting a curve
    else:
        power = -Iq(mu_CV, eG, eC, eI, eM)*mu_CV # This line is used when optimize.fminbound() calls the function
    return power

def neg_P_mu(mu_CV, eG, eI, eM, Ts = 6000, Tc = 300):
    # This function is the negative of Jv(). It is used by fminboun() to find the max of Jv()
    # I tried "return -P()" but Python it didn't like that

    eC = eG - eI

    if isinstance(mu_CV, list) or isinstance(mu_CV, np.ndarray):
        power = [-Iq(i, eG, eC, eI, eM)*i for i in mu_CV] # This line is used for plotting a curve
    else:
        power = -Iq(mu_CV, eG, eC, eI, eM)*mu_CV # This line is used when optimize.fminbound() calls the function
    return power

####################### Current curve ######################

def J(eG, mu_CV, eI, eM, Ts = 6000, Tc = 300):
    eC = eG - eI

    if isinstance(mu_CV, list) or isinstance(mu_CV, np.ndarray):
        N = len(mu_CV)
        power = np.zeros(N)
        for i in range(N):
            power[i] = Iq(mu_CV[i], eG, eC, eI, eM) # This line is used for plotting a curve
    else:
        power = Iq(mu_CV, eG, eC, eI, eM) # This line is used when optimize.fminbound() calls the function
    return power

############################################################

def optimise_muCV(eG, eI, eM, Ts = 6000, Tc = 300):
    # returns optimal mu_CV given eG and eI
    initial_mu_CV = eG*0.95
    maximum = optimize.minimize(neg_P_mu, initial_mu_CV, args = (eG, eI, eM), bounds = [(0,eG)], tol = 1e-6)
    return maximum.x[0]

############################################################

def optimise_eG(eG, eI, eM, Ts = 6000, Tc = 300):
    # Takes starting guesses of eG and returns optimal eG.

    ideal_mu_CV = optimise_muCV(eG, eI, eM, Ts, Tc)
    print(f'eG = {eG}')
    print(f'eI = {eI}')
    print(f'Ideal mu_CV: {ideal_mu_CV}')
    print('')
    return neg_P_mu(ideal_mu_CV, eG, eI, eM)

############################################################

def optimise_all(params, eM, Ts = 6000, Tc = 300):
    # Takes starting guesses of eG and eI and returns optimal eG, eI.
    eG = params[0]
    eI = params[1]

    ideal_mu_CV = optimise_muCV(eG, eI, eM, Ts, Tc)
    print(f'eG = {eG}')
    print(f'eI = {eI}')
    print(f'Ideal mu_CV: {ideal_mu_CV}')
    print('')
    return neg_P_mu(ideal_mu_CV, eG, eI, eM)

############################################################

main()
