class grid3D(object):

    """
                      #####  ######  ### ######     #     # ####### ######  #     # #       #######
                     #     # #     #  #  #     #    ##   ## #     # #     # #     # #       #
                     #       #     #  #  #     #    # # # # #     # #     # #     # #       #
                     #  #### ######   #  #     #    #  #  # #     # #     # #     # #       #####
                     #     # #   #    #  #     #    #     # #     # #     # #     # #       #
                     #     # #    #   #  #     #    #     # #     # #     # #     # #       #
                      #####  #     # ### ######     #     # ####### ######   #####  ####### #######


                    Calculates the energy grid for the cif file using a cpu, gpu or read it from a .cube file	
                    
    """

    def __init__(self, cif_file, spacing = 1.0, cutoff = 12.8):
    
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
        import numpy as np
        from pymatgen.io.cif import CifParser
        self.spacing=spacing
        self.cutoff=cutoff
        f1 = CifParser(cif_file)
        struct = f1.get_structures(primitive=False)[0]

        # Need the box lengths to ensure the distance criterion
        # Create the cell matrix to ensure this criterion

        la = struct.lattice.a
        lb = struct.lattice.b
        lc = struct.lattice.c

        alpha = struct.lattice.alpha * (np.pi/180.0)
        beta = struct.lattice.beta * (np.pi/180.0)
        gamma = struct.lattice.gamma * (np.pi/180.0)
        vol = struct.volume

        eA = [la, 0, 0]
        eB = [lb*np.cos(gamma), lb*np.sin(gamma), 0]
        eC = [lc*np.cos(beta), lc*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma), vol/(la*lb*np.sin(gamma))]

        # Find the perpendicular box lengths.
        # Those are the projections of the lattice vectors on the x, y and z axes
        # it can be shown that these lengths are equal to the inverse magnitude of the corresponding reciprocal vectors
        #  Eg . a.i = 1/|a*|

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

        struct.make_supercell([self.nx_cells, self.ny_cells, self.nz_cells])  # Structure is made into a super cell
        coord = np.array(struct.frac_coords)  # The whole thing scaled to [0,1] in all D's
        self.number_of_atoms = struct.num_sites

        # Redefine the box matrix since we made a supercell
        self.la = struct.lattice.a
        self.lb = struct.lattice.b
        self.lc = struct.lattice.c
        self.alpha = struct.lattice.alpha * (np.pi/180.0)
        self.beta = struct.lattice.beta * (np.pi/180.0)
        self.gamma = struct.lattice.gamma * (np.pi/180.0)
        self.vol = struct.volume

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

        # Read the corresponding forcefield parameters (pseudo atoms information)
        # for the atoms in the MOF

        # Get atom names in the MOF
        self.mof_atm_names = []
        for i in range(self.number_of_atoms):
                self.mof_atm_names.append(str(struct.species[i]))
        
        # def CalculateCPU():
    
