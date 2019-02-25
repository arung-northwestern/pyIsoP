class forcefields:
    """
     #######
     #        ####  #####   ####  ######    ###### # ###### #      #####   ####
     #       #    # #    # #    # #         #      # #      #      #    # #
     #####   #    # #    # #      #####     #####  # #####  #      #    #  ####
     #       #    # #####  #      #         #      # #      #      #    #      #
     #       #    # #   #  #    # #         #      # #      #      #    # #    #
     #        ####  #    #  ####  ######    #      # ###### ###### #####   ####
     * Houses the force field object for the grid calculation
     * Houses the function to read the force field from the RASPA format dir and mixing rules file
     ! Default params are for hydrogen from UFF
    """

    def __init__(self, grid_obj, force_field='UFF', sigma=2.958, epsilon=36.7):
        """ Description
            
        :type grid_obj: grid3D instance
        :param grid_obj: The grid3D object for the said cif file
    
        :type force_field: str
        :param force_field: Name of the RASPA style force field directory
    
        :type sigma: float 
        :param sigma: Lennard Jones zero potential distance
    
        :type epsilon: float 
        :param epsilon: Lennar Jones well-depth
    
        :raises:
    
        :rtype: Force field object containing all the cross interaction parameters for a given cif file after MIC
        """
        import numpy as np

        self.epsilon_array = []
        self.sigma_array = []
        self.force_field = force_field
        self.sigma = sigma
        self.epsilon = epsilon
        ff_file = open('./forcefield/'+self.force_field+'/force_field_mixing_rules.def').readlines()
        mof_atm_indices = np.zeros((grid_obj.number_of_atoms, 1))
        for i in range(len(grid_obj.mof_atm_names)):
                for line in ff_file:
                        if grid_obj.mof_atm_names[i] == line.split()[0]:
                                mof_atm_indices[i] = ff_file.index(line)

        for atm_index in range(grid_obj.number_of_atoms):

                # Lorentz Berthelot
                sig2 = float(ff_file[int(mof_atm_indices[atm_index])].split()[3])  # Get sigma
                eps2 = float(ff_file[int(mof_atm_indices[atm_index])].split()[2])  # get epsilon

                sig = (self.sigma + sig2) / 2
                eps = np.sqrt(self.epsilon * eps2)
                self.sigma_array.append(sig)
                self.epsilon_array.append(eps)

        self.sigma_array = np.array(self.sigma_array)
        self.epsilon_array = np.array(self.epsilon_array)
