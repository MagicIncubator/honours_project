import numpy as np
import matplotlib.pyplot as plt
import csv
from pandas_ods_reader import read_ods
import pandas
from pint import UnitRegistry
import warnings

U = UnitRegistry()

####################################################

def main():
    ####################################################
    ### Parameters

    h = 6.62607015e-34 * U.joule * U.second # Js
    c = 299792458       * U.meter /U.second # m/s
    q = 1.602176634e-19 * U.joule / U.eV    # J/eV
    h = h.to('eV*s') # eV s

    photon_weighted = False # if you want to calculate using I(E)/E, set to True
    plotting = True # if you want to see some plots describing each step of the calculation, set to True

    # Absorption scaling factor for the absorptions so it shows up on the plot with the spectrum data.
    scaling1 = 1 # for plotting as function of lambda
    scaling2 = -6 # for plottin as function of E

    Ts = 5780 * U.kelvin #K
    E_min = 1.9228 * U.eV # eV: is the energy values of the lowest transition
    # 1.9001388621073196
    # 1.9228
    material = 'AlGaAs'

    # if material is not an alloy, set molar_fraction = 0
    molar_fraction = 0.411 # Options for AlGaAs: 0.411, 0.419, 0.590

    ####################################################
    ### Retreiving data
    wvlgth_a, alpha_k, E, alpha_E = get_alpha_k(material, molar_fraction) # wvlght_a for "alpha"
    wvlgth_s, etr, tilt, circumsolar = get_atm1_5() # wvlght_s for "spectrum" or "sun"
    wvlgth_bb, blackbody = blackbody_distribution(Ts)

    ## Sometimes you can fined the values of alpha on a scanned pdf.
    ## If you have both k and the value of alpha lifted (copy/pasted)
    ## from a datasheet, you can compare the difference in the line below:
    #plt.figure()
    #plt.plot(wvlgth_a, np.abs(alpha*1e-3 - alpha_419)/alpha_419, label = 'Relative difference in alphas')
    #plt.show()

    #####################################################################################
    #### Plotting raw data

    if plotting:
        fig, ax1 = plt.subplots()
        plt.title('All raw data as a function of wavelength')
        plt.xlabel(r'Wavelength $\lambda$ (m)')
        ax1.set_ylabel(r"Solar spectrum intensity ($W~m^{-3}$)", color = 'blue')
        ax1.plot(wvlgth_s, circumsolar, label = "AM1.5D", color = 'blue')
        ax1.tick_params(axis ='y', labelcolor = 'blue')
        ax2 = ax1.twinx()

        ax2.set_ylabel('Absorption ($cm^{-1}$)', color = 'red')
        ax2.plot(wvlgth_a, alpha_k.to('1/cm'), label = r'$\alpha(\lambda)$ for AlGaAs, x = {}'.format(molar_fraction), color = 'red', zorder = 2)
        ax2.plot(0,1, label = "AM1.5D", color = 'blue', zorder = 1) # this was my hack for getting the label to appear
        ax2.hlines(2e4, 0, 4e-6, color = 'k', linestyle = '--', label = r'Current approximation of $\alpha(E)$')
        ax2.tick_params(axis ='y', labelcolor = 'red')
        plt.legend()

        plt.show()
    ######################################################################################
    #### jacobian transformation

    # all the same Es
    # eV, W m⁻² eV⁻¹                   meters         W m⁻³
    Es, circumsolar_E = get_jacobian(wvlgth_s, circumsolar)#*1e-9, *1e9 # eV, s⁻¹ m⁻¹ nm⁻¹
    _, tilt_E = get_jacobian(wvlgth_s, tilt)#*1e-9, *1e9
    _, etr_E = get_jacobian(wvlgth_s, etr)#*1e-9, *1e9

    if photon_weighted == True:
        # using photon-weighted
        circumsolar_E = circumsolar_E/Es
        tilt_E = tilt_E/Es
        etr_E = etr_E/Es

    if plotting:

        ## Vs E (instead of vs wavelength)
        fig, ax1 = plt.subplots()
        plt.xlabel(r'Energy (eV)')
        plt.title('All raw data as a function of energy')
        #plt.ylabel(r"$W~m^{-3}$; absorption coefficient $(cm^{-1})$") # units according to data sheet
        ax1.set_ylabel(r"Solar spectrum intensity ($W~m^{-2}~eV^{-1} $)", color = 'blue')
        ax1.plot(Es, circumsolar_E, label = "AM1.5D", color = 'blue')
        ax1.tick_params(axis ='y', labelcolor = 'blue')
        ax2 = ax1.twinx()

        ax2.set_ylabel(r'Absorption ($cm^{-1}$)', color = 'red')
        ax2.plot(E, alpha_E.to('1/cm'), label = r'$\alpha(E)$ for AlGaAs, x = {}'.format(molar_fraction), color = 'red')
        ax2.tick_params(axis ='y', labelcolor = 'red')
        ax2.plot(0,1, label = "AM1.5D", color = 'blue', zorder = 1) # this was my hack for getting the label to appear
        ax2.hlines(2e4, 0, 4.5, color = 'k', linestyle = '--', label = r'Current approximation of $\alpha(E)$')
        plt.legend()
        #ax1.legend()
        #plt.savefig('AM15D_AlGaAs_alpha_energy.png')
        #plt.show()

    ############################################
    ### integrating to check work
    integral = reimann_sum(Es, tilt_E)

    integral1 = reimann_sum(Es, tilt_E)
    integral2 = reimann_sum(Es, circumsolar_E)
    integral3 = reimann_sum(Es, etr_E)
    print('AM1.5G:', round(integral1) )#, 'W/m²')
    print('AM1.5D:', round(integral2) )#, 'W/m²')
    print('AM0   :', round(integral3) )#, 'W/m²')

    ######################################################
    ### finding the range of wavelengths over which we have both values for the solar spectrum and the absorption coefficient

    # we need to find which points overlap in the wvlth datasets
    # stripping units, becuase np doesn't play nice with pint
    e = Es.units
    c = circumsolar_E.units
    a = alpha_E.units

    Es_cut, Ea_cut, circumsolar_E_cut, alpha_E_cut = find_overlap_4(Es.magnitude[::-1], E.magnitude[::-1], circumsolar_E.magnitude[::-1], alpha_E.magnitude[::-1])

    # putting units back on
    alpha_E_cut = alpha_E_cut * a
    Es_cut = Es_cut * e
    Ea_cut = Ea_cut * e
    circumsolar_E_cut = circumsolar_E_cut * c

    if plotting:

        fig, ax1 = plt.subplots()
        plt.title('Where we have data for on both the solar spectrum and absorption')
        plt.xlabel(r'Energy (eV)')
        #plt.ylabel(r"$W~m^{-3}$; absorption coefficient $(cm^{-1})$") # units according to data sheet
        ax1.set_ylabel(r"Solar spectrum intensity ($W~m^{-2}~eV^{-1} $)", color = 'blue')
        ax1.plot(Es_cut, circumsolar_E_cut, label = "AM1.5D", color = 'blue')
        ax1.tick_params(axis ='y', labelcolor = 'blue')
        ax2 = ax1.twinx()

        ax2.set_ylabel(r'Absorption ($cm^{-1}$)', color = 'red')
        ax2.plot(Ea_cut, alpha_E_cut.to('1/cm'), label = r'$\alpha(E)$ for AlGaAs, x = {}'.format(molar_fraction), color = 'red', zorder = 2)
        ax2.plot(1,1, label = "AM1.5D", color = 'blue', zorder = 1) # this was my hack for getting the label to appear
        ax2.tick_params(axis ='y', labelcolor = 'red')
        plt.legend()
        #plt.savefig('AM15D_AlGaAs_alpha_energy_dataOnBoth.png')
        #plt.show()


    ######################################################
    ### interpolation of alpha and spectrum

    # finding the approximate value of alpha at the spectrum's wavelength points
    # Alpha is smoother than the spectrum, so it's easier to alpha at the spectrum datapoints
    # (also, the spectrum data has more points per dE).
    alpha_interp = np.interp(Es_cut, Ea_cut, alpha_E_cut)

    if plotting:
        plt.figure()
        plt.title('Comparing interpolated data with original')
        plt.plot(Ea_cut, alpha_E_cut, label = r"$\alpha_{{{}}}$".format(molar_fraction) )
        plt.plot(Es_cut, alpha_interp, label = r"Interpolated $\alpha_{{{}}}$".format(molar_fraction) )
        plt.xlabel(r'Energy (eV)')
        plt.ylabel(r"Abosroption coefficient ($eV^{-1}$)")
        plt.legend()
        #plt.close()
        #plt.show()

    ######################################################
    ### removing all values where alpha = 0 since those wavelengths can't contribute towards the generation of current
    Es_new, alpha_new, spectrum_new = remove_zeros_alpha( Es_cut, alpha_interp, circumsolar_E_cut)
    #Es_new, alpha_new, spectrum_new = Es_cut.copy(), alpha_interp.copy(), circumsolar_E_cut.copy()

    if plotting:

        fig, ax1 = plt.subplots()
        plt.title(r'All data where $\alpha \neq 0$')
        plt.xlabel(r'Energy (eV)')
        #plt.ylabel(r"$W~m^{-3}$; absorption coefficient $(cm^{-1})$") # units according to data sheet
        ax1.set_ylabel(r"Solar spectrum intensity ($W~m^{-2}~eV^{-1} $)", color = 'blue')
        ax1.plot(Es_new, spectrum_new, label = "AM1.5D", color = 'blue')
        ax1.tick_params(axis ='y', labelcolor = 'blue')
        ax2 = ax1.twinx()

        ax2.set_ylabel(r'Absorption ($cm^{-1}$)', color = 'red')
        ax2.plot(Es_new, alpha_new.to('1/cm'), label = r'$\alpha(E)$ for AlGaAs, x = {}'.format(molar_fraction), color = 'red', zorder = 2)
        ax2.plot(2,1, label = "AM1.5D", color = 'blue', zorder = 1) # this was my hack for getting the label to appear
        ax2.tick_params(axis ='y', labelcolor = 'red')
        plt.legend()
        #plt.savefig('AM15D_AlGaAs_alpha_energy_dataOnBoth.png')
        #plt.show()

    ######################################################
    ### removing all values E<E_min

    integrand_not_cut = spectrum_new*alpha_new
    E_integrand_not_cut = Es_new
    Es_new, alpha_new, spectrum_new = remove_low_E( E_min, Es_new, alpha_new, spectrum_new)
    #Es_new, alpha_new, spectrum_new = Es_new, alpha_new, spectrum_new
    integrand = spectrum_new*alpha_new


    if plotting:

        fig, ax1 = plt.subplots()
        plt.title(r'All data where $\alpha \neq 0$ & where $E>${} eV'.format(E_min.magnitude))
        plt.xlabel(r'Energy (eV)')
        ax1.set_ylabel(r"Solar spectrum intensity ($W~m^{-2}~eV^{-1} $)", color = 'blue')
        ax1.plot(Es_new, spectrum_new, label = 'cut down AM1.5D', color = 'blue')
        ax1.tick_params(axis ='y', labelcolor = 'blue')
        ax2 = ax1.twinx()

        ax2.set_ylabel(r'Absorption ($cm^{-1}$)', color = 'red')
        ax2.plot(Es_new, alpha_new.to('1/cm'), label = r'$\alpha(E)$ for AlGaAs, x = {}'.format(molar_fraction), color = 'red', zorder = 2)
        ax2.plot(2,1, label = 'cut down AM1.5D', color = 'blue', zorder = 1) # this was my hack for getting the label to appear
        ax2.tick_params(axis ='y', labelcolor = 'red')
        plt.legend()
        #plt.savefig('AM15D_AlGaAs_alpha_energy_dataOnBoth.png')
        #plt.show()


        plt.figure()
        plt.title(r'All non-zero $\alpha(E) \times I(E)$ for {} with $x={}$'.format(material, molar_fraction))
        plt.vlines(E_min.magnitude, 0, max(integrand.magnitude), color = 'k', linestyle = '--', label = 'Lowest energy required for a transition\n(aka lower bound of integration)')
        plt.plot(Es_new, integrand, c = 'b')
        plt.plot(E_integrand_not_cut, integrand_not_cut, c = 'b')
        plt.xlabel(r'Energy (eV)')
        #plt.ylabel(integrand.units)
        plt.ylabel(r'$W~eV^{-1}~m^{-3}$')
        plt.legend()
        plt.savefig('integrad.png')
        plt.grid()
        plt.show()
        #plt.close('all')

    ######################################################
    ### Final calculations -- integrals

    integrand  = alpha_new*spectrum_new
    numerator = reimann_sum(Es_new, integrand)
    denominator = reimann_sum(Es_new, spectrum_new)
    avg_absorption = (numerator/denominator).to('1/cm')
    print('\nUsing Reiman sums to calculate the integrals')
    print( 'Avg absorption(E) =', round(avg_absorption.magnitude,4), '  ', avg_absorption.units )

    ######################################################
    ##### Comparing to avg(integrad)/avg(E)
    # NOTE: don't use these results! This is for educational purposes only
    """
    avg_abs_x_spectrum = sum(integrand)/len(integrand)
    #print('avg( I(E) * alpha(E) ) =', round(avg_abs_x_spectrum.magnitude, 3), avg_abs_x_spectrum.units)
    avg_E = sum(spectrum_new)/len(spectrum_new)
    avg_absorption = avg_abs_x_spectrum/avg_E
    avg_absorption = avg_absorption.to('1/cm')
    print('\nUsing avg to calculate the integrals')
    print( 'Avg absorption =', round(avg_absorption.magnitude,4), '  ', avg_absorption.units )
    #"""

########################################################################################
########################################################################################
# DON'T EDIT PAST THIS POINT
# (unless you know what you're doing)
########################################################################################
########################################################################################

def get_atm1_5():
    base_path = "absorption&spectrum_data/AM1_5.ods"
    sheet_index = 0
    df = read_ods(base_path , sheet_index , columns=[ "Wvlgth nm", "Etr W*m-2*nm-1", "Global tilt  W*m-2*nm-1", "Direct+circumsolar W*m-2*nm-1"])

    wvlgth = df["Wvlgth nm"][1:].to_numpy().astype(float) * U.nm
    #print(wvlgth)
    etr = df["Etr W*m-2*nm-1"][1:].to_numpy().astype(float) * U.W /U.m**2 /U.nm
    tilt = df["Global tilt  W*m-2*nm-1"][1:].to_numpy().astype(float) * U.W /U.m**2 /U.nm
    circumsolar = df["Direct+circumsolar W*m-2*nm-1"][1:].to_numpy().astype(float) * U.W /U.m**2 /U.nm

    wvlgth = wvlgth.to('m')
    etr = etr.to('W/m^3')
    tilt = tilt.to('W/m^3')
    circumsolar = circumsolar.to('W/m^3')
    return wvlgth, etr, tilt, circumsolar

####################################################

def get_alpha_k(material, molar_fraction = 0):
    h = 6.62607015e-34 * U.joule * U.second # Js
    c = 299792458       * U.meter /U.second # m/s
    q = 1.602176634e-19 * U.joule / U.eV # J/eV
    h = h/q
    # https://refractiveindex.info/?shelf=other&book=AlAs-GaAs&page=Aspnes-59.0

    X = f'0{int(molar_fraction*1000)}' # molar fraction *1000 of AlGaAs, this is for getting the filename

    if material == 'AlGaAs':
        filename = f'absorption&spectrum_data/{material}/RefractiveIndexINFO_{X}.csv'
    else:
        filename = f'absorption&spectrum_data/{material}/RefractiveIndexINFO.csv'
        if material == 'AlAs' or material == 'AlN':
            warnings.warn(f'The data on {material} is not experimental. Instead it is from models and simulations. \nSource: https://refractiveindex.info/?shelf=main&book=AlAs&page=Rakic')

    df = pandas.read_csv(filename)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wvlgth = df["Wavelength, µm"].to_numpy().astype(float) * U.micrometer

    k = df["k"].to_numpy().astype(float)

    #########  alpha(lambda) #########
    # convert from um to m
    wvlgth = wvlgth.to('m')
    alpha_lam = (4*np.pi*k/wvlgth) # cm⁻¹

    #########  alpha(E) #########
    E = h*c/wvlgth
    alpha_E = 4*np.pi*k*E/(h*c)

    return wvlgth, alpha_lam, E, alpha_E

####################################################

def get_jacobian(lam, y):
    """
    Takes lam (the wavelength in meters) and a function y(lam)
    and returns the jacobian transformation into E and y(E)

    -> E = hc/lam
    -> f(E) = f(lam) * (-hc/E²)
    """

    h = 6.62607015e-34 * U.joule * U.second # Js
    c = 299792458       * U.meter /U.second # m/s
    q = 1.602176634e-19 * U.joule / U.eV    # J/eV
    h = h/q # eV s

    lam = lam.to('m') # making sure lambda is in meter if it isn't already

    E = (h*c/lam) # (eV s * (m/s) /m) = eV

    jac = (h*c/E**2) # ( (eV s * m/s) / eV² ) = m/eV
    y_E = (y*jac)
    # when applying this to the intensity vs lambda you get:
    #   W m⁻²nm⁻¹ * m/eV = W m⁻² nm⁻¹ eV⁻¹, where we convert J to eV
    #   W m⁻³ * m/eV = W m⁻² eV⁻¹, where we convert J to eV
    #
    # When applying this to the absorption coefficient vs lambda, you get:
    #   cm⁻¹ * m/J = m/cm eV⁻¹ = 100 eV⁻¹

    return E, y_E

####################################################

def blackbody_distribution(T):
    ### CURENTLY NOT USED
    """
    Takes the wavelengths (in meters) over which to plot the blackbody distribution
    Takes T, the temperature (in K) of the blackbody
    """
    #wvlgth_m = wvlgth*1e-9
    wvlgth = np.linspace(1e-9,4000e-9, 4000)*U.meter
    h = 6.62607015e-34 * U.joule * U.second # Js
    c = 299792458       * U.meter /U.second # m/s
    kB = 1.380649e-23   * U.joule / U.kelvin # J/K
    sigma = 5.670374419e-8 * U.W * U.m**(-2)*U.kelvin**(-4)

    r_earth = 6378000 * U.meter
    r_sun = 696342000 * U.meter

    fs = r_earth**2/(4*r_sun**2)

    blackbody = fs*((2*np.pi*h*c**2)/(wvlgth**5)) * 1/( np.exp(h*c/( wvlgth*kB*T ) ) -1 )# W m⁻³
    return wvlgth, blackbody

####################################################

def reimann_sum(x, f):
    if len(x) != len(f):
        raise Error(f"'x' and 'f' must be of same size. They are of size {len(x)} and {len(f)}")

    integral = 0
    for i in range(1, len(f)):
        delta_x = x[i-1] - x[i]
        integral += f[i-1] * delta_x

    return integral

####################################################

def find_overlap_4(x1, x2, y1, y2):
    """
    x1, x2: two monotonically increasing 1D arrays. find_overlap_4 returns two sub arrays
            containing the values of x1 and x2 that overlap with each other. The returned
            sub matrices x1_cut and x2_cut will be defined in such a way that
            x1_cut[0] < x2_cut[0] and x1_cut[-1] > x2_cut[-1].
    y1, y2: y1 is a function of x1 and y2 a function of x2. Both will be cut down the sam
            way as their respective domain arrays.

    Example:
    >>> x1 = [ 0.  1.  2.  3.  4.  5.  6.  7.  8.  9. 10.]
    >>> x2 = [ 6.5  7.5  8.5  9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 19.5 20.5]
    >>> y1 = [   0.    3.   12.   27.   48.   75.  108.  147.  192.  243.  300.  363.  432.  507.  588.  675.  768.  867.  972. 1083. 1200.]
    >>> y2 = [ 126.75  168.75  216.75  270.75  330.75  396.75  468.75  546.75  630.75  720.75  816.75  918.75 1026.75 1140.75 1260.75]
    >>> find_overlap(x1,x2,y1,y2)
     [ 6.  7.  8.  9. 10. 11. 12. 13. 14. 15. 16. 17. 18. 19. 20.]
     [ 6.5  7.5  8.5  9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 19.5]
     [ 108.  147.  192.  243.  300.  363.  432.  507.  588.  675.  768.  867.
       972. 1083. 1200.]
     [ 126.75  168.75  216.75  270.75  330.75  396.75  468.75  546.75  630.75
       720.75  816.75  918.75 1026.75 1140.75]
    """

    # Check that both arrays are increasing
    if x1[1] < x1[0]: # if the array is ordered in decreasing order
        x1 = x1[::-1] # reverse the arrays
        y1 = y1[::-1] # and array of y1(x1)

    if x2[1] < x2[0]: # if the array is ordered in decreasing order
        x2 = x2[::-1] # reverse the arrays
        y2 = y2[::-1] # and array of y1(x1)

    if x1[0] < x2[0]:
        i = 1
        while x1[i] < x2[0]:
            i+=1
            if i >= len(x1):
                warnings.warn(f'Arrays x1 and x2 do not overlap. The array x1 spans [{x1[0]}, {x1[-1]}] and the array x2 spans [{x2[0]}, {x2[-1]}]')
                x1_cut = []
                x2_cut = []
                y1_cut = []
                y2_cut = []
                return x1_cut, x2_cut, y1_cut, y2_cut
        i-=1
        x1_cut = x1.copy()[i:]
        x2_cut = x2.copy()
        y1_cut = y1.copy()[i:]
        y2_cut = y2.copy()

    else:
        i = 1
        while x2[i] < x1[0]:
            i+=1
            if i >= len(x2):
                warnings.warn(f'Arrays x1 and x2 do not overlap. The array x1 spans [{x1[0]}, {x1[-1]}] and the array x2 spans [{x2[0]}, {x2[-1]}]')
                x1_cut = []
                x2_cut = []
                y1_cut = []
                y2_cut = []
                return x1_cut, x2_cut, y1_cut, y2_cut
        i-=1
        x1_cut = x1.copy()
        x2_cut = x2.copy()[i:]
        y1_cut = y1.copy()
        y2_cut = y2.copy()[i:]

    # finding where they stop overlapping
    if x1_cut[-1] < x2_cut[-1]:
        i = len(x2_cut)-1
        while x2_cut[i] > x1_cut[-1]:
            i-=1
        i+=1
        x1_cut = x1_cut
        x2_cut = x2_cut[:i]
        y1_cut = y1_cut
        y2_cut = y2_cut[:i]
    else:
        i = len(x1_cut)-1
        while x1_cut[i] > x2_cut[-1]:
            i-=1
        i+=2
        x1_cut = x1_cut[:i]
        x2_cut = x2_cut
        y1_cut = y1_cut[:i]
        y2_cut = y2_cut

    return x1_cut, x2_cut, y1_cut, y2_cut

####################################################

def remove_zeros_alpha( Es, alpha, spectrum):
    alpha_new = alpha[alpha != 0]
    spectrum_new = spectrum[alpha != 0]
    Es_new = Es[alpha != 0]
    return Es_new, alpha_new, spectrum_new

####################################################

def remove_low_E(E_min, E, alpha, spectrum):
    """
    Removes: all values in E where E[i] < E_min,
           : all values alpha(E<E_min)
           : all values spectrum(E<E_min)
    Returns: the new E, alpha, spectrum
    """
    E_new = E[E>=E_min]
    alpha_new = alpha[E>=E_min]
    spectrum_new = spectrum[E>=E_min]
    return E_new, alpha_new, spectrum_new

####################################################
if __name__ == '__main__':
    main()

####################################################