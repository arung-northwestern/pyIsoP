class potentials:
    
        """
         ######
         #     #  ####  ##### ###### #    # ##### #   ##   #       ####
         #     # #    #   #   #      ##   #   #   #  #  #  #      #
         ######  #    #   #   #####  # #  #   #   # #    # #       ####
         #       #    #   #   #      #  # #   #   # ###### #           #
         #       #    #   #   #      #   ##   #   # #    # #      #    #
         #        ####    #   ###### #    #   #   # #    # ######  ####
         # TODO: Add the Feynman-Hibbs potential along with the functional form to the list
         
        """
        def __init__(self):
            pass

        def lj(self, rsq, epsilon =36.7 , sigma = 2.958): # * Defaults are for hydrogen
    
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
            self.epsilon = epsilon
            self.sigma = sigma
            return (4*self.epsilon) * (((self.sigma)**12/(rsq)**6) - ((self.sigma)**6/(rsq)**3))
        # def ljfh(self, rsq, epsilon=36.7, sigma=2.958):
        #     self.epsilon = epsilon
        #     self.sigma = sigma
        #     return (4*self.epsilon) * (((self.sigma)**12/(rsq)**6) - ((self.sigma)**6/(rsq)**3))
