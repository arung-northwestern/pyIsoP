class siter():
    """
     #####  ### ####### ####### ######
    #     #  #     #    #       #     #
    #        #     #    #       #     #
     #####   #     #    #####   ######
          #  #     #    #       #   #
    #     #  #     #    #       #    #
     #####  ###    #    ####### #     #

    * Contains the subroutine to assign molecule beads to minima inside a MOF
    
    """

    def __init__(self, grid_obj, path_to_component_pdb, ntoms_gas, spacing=1.0, model_number=0, optim_method='l-bfgs-b'):
        self.grid=grid_obj
    
        """ Description
        Initializes the siter object
        :type grid_obj: instance of the grid 3D class
        :param grid_obj: contains all the information regarding the grid
    
        :type path_to_component_pdb: str
        :param path_to_component_pdb: path to the RASPA output component file including extension    
        
        :type ntoms_gas: int
        :param ntoms_gas: number of atoms in one molecule of the guest
    
        :type spacing: float
        :param spacing: grid spacing DEFAULT: 1 A
    
        :type model_number: int
        :param model_number: configuration from the pdb file that you want to read
    
        :type optim_method: str
        :param optim_method: method for minimization from scipy.optimize.minimize
    
        :raises:
    
        :rtype: an instance of the siter class
        """
        self.comp=path_to_component_pdb
        self.nguest=ntoms_gas
        self.epsilon=98
        self.sigma=3.73
        self.spacing=spacing
        self.potential='lj'
        self.method=optim_method
        self.index=model_number

    def  molecule2minima(siter):
        
    
        """ Description
        Assigns molecule to the nearest minima for siting analyzes
        
        :type siter: instance of the grid 3D class
        :param siter:contains all the information regarding the siting analysis
    
        :raises:
    
        :rtype: assignment array, warning array (0- safe, 1-warning), shortest_dist array
        """    
        # * Interpolator with pbc
        def pot_mod(X):
            X                        = np.dot(t2.A_inv, X)  # Now we are in the unit box world
            X                        = X-np.floor(X)     # This should be floor for pbc and nint for minimum image distance
            return pot_inter(np.around(X, decimals=3))

        import numpy as np
        from scipy.optimize import minimize
        from ase.io import read, write
        import grid3D
        import forcefields
        import potentials
        import writer

        # print("Step 1: Calculating the energy grid")
        # #* Calculate the grid and find all the local minima
        # t1                         = grid3D.grid3D(siter.frame, spacing=siter.spacing) # Initialize grid
        # f1                         = forcefields.forcefields(t1, sigma=siter.sigma, epsilon=siter.epsilon) # mixing rules
        # t2                         = grid3D.grid3D.grid_calc(t1, siter.potential, f1) # Calculate grid
        #
        print("Step 1: Calculating the Local minima")
        t2                           = siter.grid
        lm                           = grid3D.grid3D.detect_local_minima(t2)  # detect minima
        min_coord                    = np.column_stack((t2.x[lm], t2.y[lm], t2.z[lm]))
        print("Step 2: Interpolating the energy grid")

        # ! Interpolate the grid over the unit box and not the actual box
        pot_inter                    = grid3D.grid3D.GridInterpolator(t2)

        print("Step 3: Reading RASPA output")

        # * Read the RASPA output
        data                         = read(siter.comp, index=siter.index)
        coord                        = data.get_positions()
        n_atoms                      = data.get_number_of_atoms()
        number_of_molecules          = int(n_atoms/siter.nguest)

        print("Step 4: Starting minimizations")
        if number_of_molecules       == 0:
            print('WARNING:No molecules in configuration '+str(siter.index))
        else:
            assign_array             = np.zeros((n_atoms, 1))  # See if we can successfully classify the beads
            warning_array            = np.zeros((n_atoms, 1))
            short_dist               = np.zeros((n_atoms, 1))
            for i in range(n_atoms):
                # * Now for each bead we will run a minimization
                res                  = minimize(pot_mod, coord[i], method='l-bfgs-b')

                # * Find the distance of the result with all the local minima, in a perfect world
                # * one these will be zero...but our world is not perfect.
                diff                 = res.x-min_coord
                diff_u               = np.dot(t2.A_inv, diff.T).T  # Otherwise the shapes won't match
                # * We can use round here because it's the difference and not the coordinate
                diff_u               = diff_u-np.round(diff_u)
                true_vec             = np.dot(t2.A, diff_u.T).T
                rsq                  = np.sum(true_vec**2, axis=1)
                assign_array[i]      = np.argsort(rsq)[0]  # * Minima at the shortest distance is saved
                short_dist[i]        = np.sort(rsq)[0] > 4.0
                if np.sort(rsq)[0] > 4.0:
                    warning_array[i] = 1

        print("Final Step(5): Printing the assignments and warnings")
        # symbols=np.reshape(np.array(['C']*n_atoms).T,(n_atoms,1))
        #out=np.hstack((coord, assign_array,warning_array))
        #np.savetxt(str(siter.index)+'_assign.xyz',out)
        return assign_array, warning_array, short_dist
