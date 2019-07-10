class predictors:
    
    #*######  ######  ####### ######  ###  #####  ####### ####### ######
    #*#     # #     # #       #     #  #  #     #    #    #     # #     #
    #*#     # #     # #       #     #  #  #          #    #     # #     #
    #*######  ######  #####   #     #  #  #          #    #     # ######
    #*#       #   #   #       #     #  #  #          #    #     # #   #
    #*#       #    #  #       #     #  #  #     #    #    #     # #    #
    #*#       #     # ####### ######  ###  #####     #    ####### #     #
    """
    * Takes in a set of pressures, void fraction, energy histogram and the first-shell coordination (n1) and predicts the adsorption isotherm
    """
    import numpy as np
    def __init__(self):
        pass

    def predict_isotherm( T, pressures, vf, hist_obj, n1, epsilon=36.7, MA=2.0, theta=0.0):
        """ 
        
        Calculates the adsorption isotherm using the fundamental equation of adsorption with the pressures, energy histogram, the void fraction
        of the material and the first shell coordination numbers of hydrogen in the material at the pressures of interest, as the inputs.
        The first shell coordination number n1 maybe calculated from the coordination module using the machine learning model, with training data
        provided only for 77 K and 160 K. Training data has to be generated using GCMC for other temperatures.
        
        :type T: float  
        :param T: Adsorption temperature

        :type pressures: numpy array of floats
        :param pressures: Adsorption pressures in bar
    
        :type vf: float 
        :param vf: void fraction of the material, needed for the saturation loading
    
        :type histogram: instance of the histo class
        :param histogram: The histogram object containing all the information about the energy histogram
    
        :type n1: Array of floats
        :param n1: First shell coordination numbers at the different pressures of interest.
    
        :type theta: Float
        :param theta: The guest-guest scaling coefficient in the Temkin approximation to BET. DEFAULT: Theta=0 ==> Langmuir form.
    
        :type epsilon: Float
        :param epsilon: The Lennard Jones well-depth (epsilon/kB) for the guest, DEFAULT: 36.7 [K] for hydrogen form UFF 
    

        :type MA: Float
        :param MA: Molecular weight of hydrogen
    
    
        :raises: None
    
        :rtype: An array of adsorption loading in g/L
        """
        import scipy.integrate as intg
        import numpy as np
        R=0.0831
        saturation_loading = 70.8 * vf
        Ered               = hist_obj.E # This E should be dimensionless
        RhoE               = hist_obj.RhoE # Should be consistent with dimesionless energy
        loading_predicted  = []
        for p_index in range(len(pressures)):
            p              = pressures[p_index]
            n1_p           = n1[p_index]
            kG             = np.exp(-(-n1_p*epsilon) / (2*T))
            kE             = np.exp(-Ered) / (R * T) / saturation_loading * MA
            #kE            = np.exp(-Ered) / (R * T) / 70.8 * MA
            k_Lang         = kG * kE
            f              = k_Lang * p / (1.0 + k_Lang * p)
            func           = saturation_loading * RhoE * (f + theta * (f**2) * (f - 1.0))
            #func          = 70.8 * Rho_E3 * (f + theta * (f**2) * (f - 1.0))
            Intg_Value     = intg.simps(func, Ered)
            loading_predicted.append(Intg_Value)
        loading_predicted  = np.reshape(loading_predicted, (len(pressures), 1))
        return loading_predicted



