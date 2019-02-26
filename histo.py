class histo:
    """
    ####### #     # ####### ######   #####  #     #    #     # ###  #####  ####### #######  #####  ######     #    #     #
    #       ##    # #       #     # #     #  #   #     #     #  #  #     #    #    #     # #     # #     #   # #   ##   ##
    #       # #   # #       #     # #         # #      #     #  #  #          #    #     # #       #     #  #   #  # # # #
    #####   #  #  # #####   ######  #  ####    #       #######  #   #####     #    #     # #  #### ######  #     # #  #  #
    #       #   # # #       #   #   #     #    #       #     #  #        #    #    #     # #     # #   #   ####### #     #
    #       #    ## #       #    #  #     #    #       #     #  #  #     #    #    #     # #     # #    #  #     # #     #
    ####### #     # ####### #     #  #####     #       #     # ###  #####     #    #######  #####  #     # #     # #     #
    * 
    TODO: Add grid object to histogram  object function
    TODO: Add text file to histogram object function
    TODO: Add RASPA grid to histogram  object function
    TODO: Add cube file to histogram object function
    ! Here the histogram is in the units of kT or RT maybe do it in [K] or kJ/mol instead
    ! Make sure to modify the predictor to account for the units of the energy in the histogram
    ! Make sure the predictor accounts for type of the histogram argument
    """
    import numpy as np

    def __init__(self, nbins=50, E_max=10, normed_flag=True):
        """ Description
        Initializes the histogram object to store the energy histogram, unit of [K] for the energy.

        :type number: int
        :param number: number of bins
    
        :type E: array of floats
        :param E: Energy of each bin in the histogram-- bincenters
    
        :type RhoE: array of floats
        :param RhoE: Normalized energy density, (area under the curve =1)

        :type E_max: float
        :param E: Maximum energy where the histogram is cut off. 
        * DEFAULT: 10 [K], doesn't change the isotherm from prior testing
    

        :raises:
    
        :rtype: Instance of the histogram class
        """

        self.nbins                                                    = number
        self.E                                                        = []
        self.RhoE                                                     = []
        self.E_max                                                    = E_max
        self.normed_flag                                              = normed_flag

    def grid2histo(grid_obj, hist_obj, T):
    
        """ Description
        :type grid_obj: An instance of the grid3D class
        :param grid_obj: Contains energy info on a 3D grid
    
        :type hist_obj: An instance of the histo class
        :param hist_obj: Will be modified and overwritten with energy and density of sites info
        ! Remove this preferably
        :type T: float   
        :param T: Temperature
    
        :raises:
    
        :rtype: Modified histogram object with E and RhoE.
        
        """    
        e_vals                                                        = np.reshape(grid_obj.pot_repeat, (grid_obj.N_grid_total, 1), order='C')
        bins1                                                         = np.linspace(min(e_vals/T)-0.5, hist_obj.E_max, hist_obj.nbins+1)
        RhoE, binedges1                                               = np.histogram(e_vals / T, bins=bins1, normed=hist_obj.normed_flag)
        bincenters                                                    = 0.5 * (binedges1[1:] + binedges1[:-1])  # Bincenters
        hist_obj.E                                                    = bincenters
        return hist_obj
