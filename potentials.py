class potentials:

        """
         ######
         #     #  ####  ##### ###### #    # ##### #   ##   #       ####
         #     # #    #   #   #      ##   #   #   #  #  #  #      #
         ######  #    #   #   #####  # #  #   #   # #    # #       ####
         #       #    #   #   #      #  # #   #   # ###### #           #
         #       #    #   #   #      #   ##   #   # #    # #      #    #
         #        ####    #   ###### #    #   #   # #    # ######  ####
      
         
        """
        from numba import vectorize
        def __init__(self):
            pass
            
        # * The Lennard-Jones Potential
        @vectorize(['float64(float64, float64, float64)'], target='cpu')
        def lj( rsq, epsilon, sigma): 
            """ Description

            :type rsq: float
            :param rsq: squared distance between atoms

            :type epsilon: float 
            :param epsilon: Depth of the potential well in [K]

            :type sigma: float
            :param sigma: Zero-energy distance in angstrom

            :raises: 
            ValueError: rsq cannot be zero

            :rtype: The Lennard-Jones potential value in 'float'.
            """
            return (4*epsilon) * (((sigma)**12/(rsq)**6) - ((sigma)**6/(rsq)**3))

        # * Lennard Jones with the Feynman-Hibbs correction    
        @vectorize(['float64(float64, float64, float64, float64, float64)'], target='cpu')
        def ljfh( rsq, epsilon, sigma, rmass, T):
            """ Description

            :type rsq: float
            :param rsq: squared distance between atoms

            :type epsilon: float 
            :param epsilon: Depth of the potential well in [K]

            :type sigma: float
            :param sigma: Zero-energy distance in angstrom
            
            :type rmass: float
            :param rmass: Reduced mass of the pair system m1*m2/(m1+m2)
            
            :type T: float
            :param T: Temperature in [K]
            
            :raises: 
            ValueError: rsq cannot be zero
            TypeError: Provide T and rmass for Feynman Hibbs correction

            :rtype: The Lennard-Jones potential value in 'float'.
            """

  
            # * Let's convert everything into SI unites to evaluate the Feynman Hibbs constant clearly: hcut**2/[(r**2)(24*pi*m*k_B*T)]
            hcut = 1.0545718E-34  # * Planck's constant / 2*pi in Js or kg m2 s-1
            angs = 1E-10          # * angstrom in m
            amu = 1.66054E-27    # * amu in kg
            k_B = 1.38065852E-27  # * Boltzmann constant in J/K or kg m2 s-2 K-1

            const = hcut**2 / ((rsq*angs**2)*(24*rmass*amu*k_B*T)) # * Note  that this constant is dimensionless

            return (4*epsilon) * (((sigma)**12/(rsq)**6) - ((sigma)**6/(rsq)**3)) + const*(4*epsilon) * (132*((sigma)**12/(rsq)**6) - 30*((sigma)**6/(rsq)**3))
