class grid3D(object):

    """

      #####  ######  ### ######     #     # ####### ######  #     # #       #######
    #     # #     #  #  #     #    ##   ## #     # #     # #     # #       #
    #       #     #  #  #     #    # # # # #     # #     # #     # #       #
    #  #### ######   #  #     #    #  #  # #     # #     # #     # #       #####
    #     # #   #    #  #     #    #     # #     # #     # #     # #       #
    #     # #    #   #  #     #    #     # #     # #     # #     # #       #
      #####  #     # ### ######     #     # ####### ######   #####  ####### #######


    Calculates the energy grid for the cif file using a cpu or gpu 
    ! Does not include the electrostatics
    TODO: Add the GPU based grid calculation..? or should we just do that as a decorator for the potential.
    """
    import numpy as np

    def __init__(self, cif_file, spacing=1.0, cutoff=12.8):
        """ Description
        Function to initiate the grid3D class object.

        :type cif_file: str
        :param cif_file: Path the cif file for which we want to calculate the grid
    
        :type spacing: float
        :param spacing: Grid spacing in angstrom, same in all directions, the default is 1.0
    
        :type cutoff: float
        :param cutoff: The LJ cut off distance in angstrom, the default is 12.8
    
        :raises:
        ValueError: The cut-off and the spacing cannot be zero.
    
        :rtype: Returns an instance of the class grid3D
        """

        from pymatgen.io.cif import CifParser
        self.file                              = cif_file
        self.spacing                           = spacing
        self.cutoff                            = cutoff
        f1                                     = CifParser(cif_file)
        struct                                 = f1.get_structures(primitive=False)[0]

        # Need the box lengths to ensure the distance criterion
        # Create the cell matrix to ensure this criterion

        la                                     = struct.lattice.a
        lb                                     = struct.lattice.b
        lc                                     = struct.lattice.c

        alpha                                  = struct.lattice.alpha * (np.pi/180.0)
        beta                                   = struct.lattice.beta * (np.pi/180.0)
        gamma                                  = struct.lattice.gamma * (np.pi/180.0)
        vol                                    = struct.volume

        eA                                     = [la, 0, 0]
        eB                                     = [lb*np.cos(gamma), lb*np.sin(gamma), 0]
        eC                                     = [lc*np.cos(beta), lc*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), vol/(la*lb*np.sin(gamma))]

        # Find the perpendicular box lengths.
        # Those are the projections of the lattice vectors on the x, y and z axes
        # it can be shown that these lengths are equal to the inverse magnitude of the corresponding reciprocal vectors
        #  Eg . a.i                            = 1/|a*|

        self.lx_unit                           = vol / np.linalg.norm(np.cross(eB, eC))
        self.ly_unit                           = vol / np.linalg.norm(np.cross(eC, eA))
        self.lz_unit                           = vol / np.linalg.norm(np.cross(eA, eB))

        # Define grid points in a unit box (Non adaptive grid)
        # Number of grid points
        self.nx                                = int(la/self.spacing)
        self.ny                                = int(lb/self.spacing)
        self.nz                                = int(lc/self.spacing)

        self.nx_cells                          = int(np.ceil(2.0*self.cutoff/self.lx_unit))    # magic formula
        self.ny_cells                          = int(np.ceil(2.0*self.cutoff/self.ly_unit))
        self.nz_cells                          = int(np.ceil(2.0*self.cutoff/self.lz_unit))

        struct.make_supercell([self.nx_cells, self.ny_cells, self.nz_cells])  # Structure is made into a super cell
        self.coord                             = np.array(struct.frac_coords)  # The whole thing scaled to [0,1] in all D's
        self.number_of_atoms                   = struct.num_sites

        # Redefine the box matrix since we made a supercell
        self.la                                = struct.lattice.a
        self.lb                                = struct.lattice.b
        self.lc                                = struct.lattice.c
        self.alpha                             = struct.lattice.alpha * (np.pi/180.0)
        self.beta                              = struct.lattice.beta * (np.pi/180.0)
        self.gamma                             = struct.lattice.gamma * (np.pi/180.0)
        self.vol                               = struct.volume

        self.A                                 = [[self.la, 0, 0], [self.lb*np.cos(gamma), self.lb*np.sin(gamma), 0], [self.lc*np.cos(beta), self.lc *
                                                                                                     (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), self.vol/(self.la*self.lb*np.sin(gamma))]]
        self.A                                 = np.array(self.A)
        self.A                                 = self.A.T
        self.A_inv                             = np.linalg.inv(self.A)

        # Intially the grids are defined only on the unit cell, which is only a tiny part of the unit box
        # The unit box corresponds to the entire super cell
        self.x_grid                            = np.linspace(0, 1.0/self.nx_cells, self.nx)
        self.y_grid                            = np.linspace(0, 1.0/self.ny_cells, self.ny)
        self.z_grid                            = np.linspace(0, 1.0/self.nz_cells, self.nz)

        # Shift every thing to the center of the box for the round commmand to work properly
        self.x_grid, self.y_grid, self. z_grid = self.x_grid - 0.5, self._grid - 0.5, self.z_grid - 0.5
        self.coord                             = self.coord - 0.5

        # Read the corresponding  atoms information
        # for the atoms in the MOF

        # Get atom names in the MOF
        self.mof_atm_names                     = []
        for i in range(self.number_of_atoms):
                self.mof_atm_names.append(str(struct.species[i]))

        self.pot                               = np.zeros((self.nx, self.ny, self.nz))
        self.pot_repeat                        = np.tile(self.pot, (self.nx_cells, self.ny_cells, self.nz_cells))
        # def CalculateCPU():

    def grid_calc(grid_obj, potential_name, ff_obj, rmass=None, T=None):

          """ Description
          :type grid_obj: instance of the grid object from the grid3D module
          :param grid_obj: contains al the details for the energy grid calculation

          :type potential_name: str
          :param potential_name: Name of the potential 'lj' or 'ljfh'

          :type ff_obj: Instance of the force field object from the forcefields module
          :param ff_obj: Contains all the mixing rules for a give functional form

          :raises:
          TypeError: If T and reduced mass for the Feynman-Hibbs is missing
          NameError: If unknown potential name is called for

          :rtype: Returns the 3D energy grid
          """

          import potentials
          p = potentials()
          potential_form = getattr(p,potential_name) # * Used the string input to choose which function to call. This is cool
          for k in range(grid_obj.nz):
            for j in range(grid_obj.ny):
              for i in range(grid_obj.nx):
                grid_point                     = np.array([grid_obj.x_grid[i], grid_obj.y_grid[j], grid_obj.z_grid[k]])
                diff_vec                       = grid_obj.coord-grid_point
                diff_vec                       = diff_vec-np.round(diff_vec)
                pbc_vec                        = np.dot(A, diff_vec.T).T
                rsq                            = np.sum(pbc_vec**2, axis=1)
                if potential_name=='lj':
                  E                              = np.sum(potential_form(ff_obj.eps_array, ff_obj.sig_array, rsq))
                elif potential_name=='ljfh' and T is not None and rmass is not None:
                  E = np.sum(potential_form(ff_obj.eps_array, ff_obj.sig_array, rsq, rmass, T))
                elif potential_name=='ljfh' and ((T is None) or (rmass is None)):
                      raise TypeError('Missing T and rmass for Feynman-Hibbs correction')
                else:
                      raise NameError('Unknown potential: Add yours to potentials.py if need be')

                grid_obj.pot[i][j][k]          = E
        

          grid_obj.pot_repeat = np.tile((grid_obj.nx_cells, grid_obj.ny_cells, grid_obj.nz_cells,),grid_obj.pot)
