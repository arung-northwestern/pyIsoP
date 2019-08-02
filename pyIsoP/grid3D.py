class grid3D(object):
    ##*    #####  ######  ### ######     #     # ####### ######  #     # #       #######
    ##*  #     # #     #  #  #     #    ##   ## #     # #     # #     # #       #
    ##*  #       #     #  #  #     #    # # # # #     # #     # #     # #       #
    ##*  #  #### ######   #  #     #    #  #  # #     # #     # #     # #       #####
    ##*  #     # #   #    #  #     #    #     # #     # #     # #     # #       #
    ##*  #     # #    #   #  #     #    #     # #     # #     # #     # #       #
    ##*    #####  #     # ### ######     #     # ####### ######   #####  ####### #######

    """
      Calculates the energy grid for the cif file using a cpu or gpu. 
      ! Does not include the electrostatics

    """

    # * Initialize the grid object
    def __init__(self, cif_file, spacing=1.0, cutoff=12.8, temperature=77.0, index=None):
        """ 
        Function to initiate the grid3D class object.

        :type cif_file: str
        :param cif_file: Path the cif file for which we want to calculate the grid
    
        :type spacing: float
        :param spacing: Grid spacing in angstrom, same in all directions, the default is 1.0
    
        :type cutoff: float
        :param cutoff: The LJ cut off distance in angstrom, the default is 12.8

        :type temperature: float
        :param temperature: Temperature at which grid is calculated. Grid has T dependence for Feynman-Hibbs and for polyatomic probe.
    

        :type index: int
        :param index: model number if you're using .pdb file
        
        :raises: ValueError
    
        :rtype: Returns an instance of the class grid3D
        
        """
        import numpy as np
        # from pymatgen.io.cif import CifParser
        from ase.io import read, write
        self.file = cif_file
        self.spacing = spacing
        self.cutoff = cutoff
        #f1                                     = CifParser(cif_file)
        frame = read(cif_file, index=None)
        #struct                                 = f1.get_structures(primitive=False)[0]

        # Need the box lengths to ensure the distance criterion
        # Create the cell matrix to ensure this criterion

        # la                                     = struct.lattice.a
        # lb                                     = struct.lattice.b
        # lc                                     = struct.lattice.c

        # alpha                                  = struct.lattice.alpha * (np.pi/180.0)
        # beta                                   = struct.lattice.beta * (np.pi/180.0)
        # gamma                                  = struct.lattice.gamma * (np.pi/180.0)
        # vol                                    = struct.volume
        # Need the box lengths to ensure the distance criterion
        # Create the cell matrix to ensure this criterion

        la = frame.get_cell_lengths_and_angles()[0]
        lb = frame.get_cell_lengths_and_angles()[1]
        lc = frame.get_cell_lengths_and_angles()[2]
        alpha = frame.get_cell_lengths_and_angles()[3] * (np.pi/180.0)
        beta = frame.get_cell_lengths_and_angles()[4] * (np.pi/180.0)
        gamma = frame.get_cell_lengths_and_angles()[5] * (np.pi/180.0)
        vol = frame.get_volume()
        eA = [la, 0, 0]
        eB = [lb*np.cos(gamma), lb*np.sin(gamma), 0]
        eC = [lc*np.cos(beta), lc*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), vol/(la*lb*np.sin(gamma))]

        # Find the perpendicular box lengths.
        # Those are the projections of the lattice vectors on the x, y and z axes
        # it can be shown that these lengths are equal to the inverse magnitude of the corresponding reciprocal vectors
        #  Eg . a.i                            = 1/|a*|

        self.lx_unit = vol / np.linalg.norm(np.cross(eB, eC))
        self.ly_unit = vol / np.linalg.norm(np.cross(eC, eA))
        self.lz_unit = vol / np.linalg.norm(np.cross(eA, eB))

        # Define grid points in a unit box (Non adaptive grid)
        # Number of grid points
        self.nx = int(la/self.spacing)
        self.ny = int(lb/self.spacing)
        self.nz = int(lc/self.spacing)

        self.nx_cells = int(np.ceil(2.0*self.cutoff/self.lx_unit))    # magic formula
        self.ny_cells = int(np.ceil(2.0*self.cutoff/self.ly_unit))
        self.nz_cells = int(np.ceil(2.0*self.cutoff/self.lz_unit))

        frame_repeat = frame.repeat([self.nx_cells, self.ny_cells, self.nz_cells])
        #struct.make_supercell([self.nx_cells, self.ny_cells, self.nz_cells])  # Structure is made into a super cell

        #self.coord                             = np.array(struct.frac_coords)  # The whole thing scaled to [0,1] in all D's
        self.coord = frame_repeat.get_scaled_positions()  # The whole thing scaled to [0,1] in all D's

        #self.number_of_atoms                   = struct.num_sites
        self.number_of_atoms = frame_repeat.get_number_of_atoms()

        # Redefine the box matrix since we made a supercell
        # self.la                                = struct.lattice.a
        # self.lb                                = struct.lattice.b
        # self.lc                                = struct.lattice.c
        # self.alpha                             = struct.lattice.alpha * (np.pi/180.0)
        # self.beta                              = struct.lattice.beta * (np.pi/180.0)
        # self.gamma                             = struct.lattice.gamma * (np.pi/180.0)
        # self.vol                               = struct.volume

        self.la = frame_repeat.get_cell_lengths_and_angles()[0]
        self.lb = frame_repeat.get_cell_lengths_and_angles()[1]
        self.lc = frame_repeat.get_cell_lengths_and_angles()[2]
        self.alpha = frame_repeat.get_cell_lengths_and_angles()[3] * (np.pi/180.0)
        self.beta = frame_repeat.get_cell_lengths_and_angles()[4] * (np.pi/180.0)
        self.gamma = frame_repeat.get_cell_lengths_and_angles()[5] * (np.pi/180.0)
        self.vol = frame_repeat.get_volume()

        self.A = [[self.la, 0, 0], [self.lb*np.cos(gamma), self.lb*np.sin(gamma), 0], [self.lc*np.cos(beta), self.lc *
                                                                                       (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), self.vol/(self.la*self.lb*np.sin(gamma))]]
        self.A = np.array(self.A)
        self.A = self.A.T
        self.A_inv = np.linalg.inv(self.A)

        # Intially the grids are defined only on the unit cell, which is only a tiny part of the unit box
        # The unit box corresponds to the entire super cell
        self.x_grid = np.linspace(0, 1.0/self.nx_cells, self.nx)
        self.y_grid = np.linspace(0, 1.0/self.ny_cells, self.ny)
        self.z_grid = np.linspace(0, 1.0/self.nz_cells, self.nz)

        # Shift every thing to the center of the box for the round commmand to work properly
        self.x_grid, self.y_grid, self. z_grid = self.x_grid - 0.5, self.y_grid - 0.5, self.z_grid - 0.5
        self.coord = self.coord - 0.5

        # Read the corresponding  atoms information
        # for the atoms in the MOF

        # Get atom names in the MOF
        self.mof_atm_names = []
        for i in range(self.number_of_atoms):
                # self.mof_atm_names.append(str(struct.species[i]))
                self.mof_atm_names.append(frame_repeat.get_chemical_symbols()[i])

        self.pot = np.zeros((self.nx, self.ny, self.nz))
        self.pot_repeat = np.tile(self.pot, (self.nx_cells, self.ny_cells, self.nz_cells))
        # * Total points
        self.nx_total = int(self.nx*self.nx_cells)
        self.ny_total = int(self.ny*self.ny_cells)
        self.nz_total = int(self.nz*self.nz_cells)

        self.N_grid_total = self.nx_total * self.ny_total * self.nz_total

        # Create the xgrid ygrid and zgrid for future use
        dx, dy, dz = 1.0/self.nx_total, 1.0/self.ny_total, 1.0/self.nz_total
        # self.X = np.arange(0, 1 + 0.1*dx, dx, dtype='float64')
        # self.Y = np.arange(0, 1 + 0.1*dy, dy, dtype='float64')
        # self.Z = np.arange(0, 1 + 0.1*dz, dz, dtype='float64')
        self.X = np.linspace(0, 1, self.nx_total)
        self.Y = np.linspace(0, 1, self.ny_total)
        self.Z = np.linspace(0, 1, self.nz_total)
        self.x = np.zeros((self.nx_total, self.ny_total, self.nz_total))
        self.y = np.zeros((self.nx_total, self.ny_total, self.nz_total))
        self.z = np.zeros((self.nx_total, self.ny_total, self.nz_total))
        for k in range(self.nz_total):
                for j in range(self.ny_total):
                        for i in range(self.nx_total):
                                self.x[i, j, k] = self.X[i]
                                self.y[i, j, k] = self.Y[j]
                                self.z[i, j, k] = self.Z[k]
        for k in range(self.nz_total):
                for j in range(self.ny_total):
                        for i in range(self.nx_total):
                                [self.x[i, j, k], self.y[i, j, k], self.z[i, j, k]] = np.dot(
                                    self.A, [self.x[i, j, k], self.y[i, j, k], self.z[i, j, k]])
        self.ase = frame_repeat
        self.pot_sphere = []
        self.pot_total = []
        self.Temperature = temperature

    # * Calculate the energy grid
    def grid_calc(grid_obj, potential_name, ff_obj, rmass=None, T=None):
          """ 
          Calculates the energy grid using 'numba' , vectorized over the cpu

          :type grid_obj: instance of the grid object from the grid3D module.
          :param grid_obj: contains al the details for the energy grid calculation.

          :type potential_name: str
          :param potential_name: Name of the potential 'lj' or 'ljfh'

          :type ff_obj: Instance of the force field object from the forcefields module
          :param ff_obj: Contains all the mixing rules for a give functional form

          :raises: TypeError, NameError
          

          :rtype: Returns the 3D energy grid
          """

       
          import numpy as np
          import pyIsoP.potentials as potentials
          from tqdm import tqdm
          p = potentials.potentials()
          # * Used the string input to choose which function to call. This is cool
          #print(potential_name)
          potential_form = getattr(p, str(potential_name))
          for k in tqdm(range(grid_obj.nz), desc='Calculating grid: Current Z position'):
            for j in range(grid_obj.ny):
              for i in range(grid_obj.nx):
                grid_point = np.array([grid_obj.x_grid[i], grid_obj.y_grid[j], grid_obj.z_grid[k]])
                diff_vec = grid_obj.coord-grid_point
                diff_vec = diff_vec-np.round(diff_vec)
                pbc_vec = np.dot(grid_obj.A, diff_vec.T).T
                rsq = np.sum(pbc_vec**2, axis=1)
                if potential_name == 'lj':
                  E = np.sum(potential_form(rsq, ff_obj.epsilon_array, ff_obj.sigma_array))
                elif potential_name == 'ljfh' and T is not None and rmass is not None:
                  E = np.sum(potential_form(rsq, ff_obj.epsilon_array, ff_obj.sigma_array, rmass, T))
                elif potential_name == 'ljfh' and ((T is None) or (rmass is None)):
                      raise TypeError('Missing T and rmass for Feynman-Hibbs correction')
                else:
                      raise NameError('Unknown potential: Add yours to potentials.py if need be')

                grid_obj.pot[i][j][k] = E

          grid_obj.pot_repeat = np.tile(grid_obj.pot, (grid_obj.nx_cells, grid_obj.ny_cells, grid_obj.nz_cells,))
          return grid_obj

    # * Interpolate the energy grid using scipy's RegularGridInterpolator
    def GridInterpolator(grid_obj):
      """ 
      Interpolates the energy grid based on scipy's RegularGridInterpolator with default settings.

      :type grid_obj: instance of the grid3D class
      :param grid_obj: contains all the information regarding the energy grid used in isotherm prediction
    
      :raises: ValueError
    
      :rtype: An interpolator function of the type RegularGridInterpolator from scipy.interpolate
      """
      from scipy.interpolate import RegularGridInterpolator
      # ! The interpolation is over a unit box, to use it on x, y, z transform by A_inv before using the interp
      pot_interp = RegularGridInterpolator((grid_obj.X, grid_obj.Y, grid_obj.Z), grid_obj.pot_repeat)
      return pot_interp

    # * Most minimization algorithms don't have PBC which is essential for us, add it in the function

    # * Detect the local minima in the energy grid, this can be used for molecular siting
    def detect_local_minima(grid_obj):

          # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
          """
          Takes an array and detects the troughs using the local maximum filter.
          Returns a boolean mask of the troughs (i.e. 1 when
          the pixel's value is the neighborhood maximum, 0 otherwise)
          This can be very useful in molecular siting later on.
          """
          import scipy.ndimage.filters as filters
          import scipy.ndimage.morphology as morphology
          import numpy as np

          arr = grid_obj.pot_repeat

          # define an connected neighborhood
          # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
          neighborhood = morphology.generate_binary_structure(len(arr.shape), 2)
          # apply the local minimum filter; all locations of minimum value
          # in their neighborhood are set to 1
          # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
          local_min = (filters.minimum_filter(arr, footprint=neighborhood) == arr)
          # local_min is a mask that contains the peaks we are
          # looking for, but also the background.
          # In order to isolate the peaks we must remove the background from the mask.
          #
          # we create the mask of the background
          background = (arr == 0)
          #
          # a little technicality: we must erode the background in order to
          # successfully subtract it from local_min, otherwise a line will
          # appear along the background border (artifact of the local minimum filter)
          # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
          eroded_background = morphology.binary_erosion(
              background, structure=neighborhood, border_value=1)
          #
          # we obtain the final mask, containing only peaks,
          # by removing the background from the local_min mask
          detected_minima = local_min - eroded_background
          return np.where(detected_minima)

    def find_minima(grid_obj):
      """ 
      Finds the locations and energy of the local minima from a pretabulated energy grid object using the detect_local_minima function.
      
      :type grid_obj: instance of the grid3D class
      :param grid_obj: contains all the information regarding the energy grid 

      :rtype: coordinates of the minima (array of floats), energy of the minima (array of floats)
      """
      import numpy as np
      import pyIsoPgrid3D as grid3D

      lm = grid3D.grid3D.detect_local_minima(grid_obj)
      out_coord = np.column_stack((grid_obj.x[lm], grid_obj.y[lm], grid_obj.z[lm]))
      pot_interp = grid3D.grid3D.GridInterpolator(grid_obj)
      en = []
      for i in range(len(out_coord)):
        point = np.dot(grid_obj.A_inv, out_coord[i])
        en.append(pot_interp(np.around(point, decimals=3)))
      #en_neg=en[en<0]
      #min_neg=out_coord[en<0]

      return np.array(out_coord), np.array(en)
