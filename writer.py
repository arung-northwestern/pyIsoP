class writer:
    """
    #     # ######  ### ####### ####### ######
    #  #  # #     #  #     #    #       #     #
    #  #  # #     #  #     #    #       #     #
    #  #  # ######   #     #    #####   ######
    #  #  # #   #    #     #    #       #   #
    #  #  # #    #   #     #    #       #    #
    ## ##  #     # ###    #    ####### #     #
    
    * Contains all the functions required to write out the details of the grid, vtk, coordinates in xyz, isotherm data and whatever
    * Most of these functions operate on the grid3D object after computation.
    ! Check citation for xyz.py file 

    """
    import os
    import imp
    import numpy as np

    xyz_mod                                                          = imp.load_source('xyz', 'xyz.py')  # A single file to make xyz file writing easy

    def __init__(self, *args, **kwargs):
        pass

    # * Details of the energy grid
    def write_details(grid_obj, path_to_file):
        """ Description
        Reads in the grid3D object after calculation and writes out the important details into a named file
        :type grid_obj: An instance of the class grid3D
        :param grid_obj: grid3D object after calculation

        :type path_to_folder: str
        :param path_to_folder: Path to the folder where you want the file to be written to.

        :raises:

        :rtype:
        """

        f1                                                           = open(path_to_file, 'w')
        f1.write(str(grid_obj.nx_total)+'\t'+str(grid_obj.ny_total)+'\t'+str(grid_obj.nz_total)+'\n')
        f1.write(str(grid_obj.nx_cells)+'\t'+str(grid_obj.ny_cells)+'\t'+str(grid_obj.nz_cells)+'\n')
        f1.write(str(grid_obj.la)+'\t'+str(grid_obj.lb)+'\t'+str(grid_obj.lc)+'\n')
        f1.write(str(grid_obj.alpha)+'\t'+str(grid_obj.beta)+'\t'+str(grid_obj.gamma)+'\n')
        f1.close()

    # * Framework coordinates in xyz
    def write_xyz(grid_obj, path_to_file):
        """ Description
        Reads in the grid3D object after calculation and writes out the important details into an xyz file for visualization
        :type grid_obj: An instance of the class grid3D
        :param grid_obj: grid3D object after calculation

        :type path_to_folder: str
        :param path_to_folder: Path to the folder where you want the file to be written to.

        :raises:

        :rtype:
        """

        import numpy as np
        f1                                                           = open(path_to_file, 'w')
        out_coord                                                    = np.zeros((np.shape(grid_obj.coord)))
        out_coord                                                    = np.dot(A, grid_obj.coord.T).T
        xyz_mod.write_xyz(f1, out_coord, title=stripped+'.xyz', atomtypes=grid_obj.mof_atm_names)
        f1.close()

    # * Energy values as an array
    def write_energy_values(grid_obj, path_to_file):
        import numpy as np

    """ Description
    :type grid_obj: instance of the grid3D class
    :param grid_obj: Contains all the info regarding the energy grid

    :type path_to_file: str
    :param path_to_file: path to the output file including full file name and extension.

    :raises:

    :rtype: Just writes the file.
    """
        f1                                                           = open(path_to_file, 'wb')  #
        e_vals                                                       = np.reshape(grid_obj.pot_repeat, (grid_obj.N_grid_total, 1), order='C')  # grid to a vector of values
        np.savetxt(f1, e_vals)
        f1.close()

    # * Energy values in x, y, z, E
    def write_xyzE(grid_obj, path_to_file):
        """ Description
        :type grid_obj: instance of the grid3D class
        :param grid_obj: Contains all the info regarding the energy grid

        :type path_to_file: str
        :param path_to_file: path to the output file including full file name and extension.

        :raises:

        :rtype: Just writes the file.
        """

        import numpy as np

        e_vals                                                       = np.reshape(grid_obj.pot_repeat, (grid_obj.N_grid_total, 1), order='C')
        grid_obj.x_total, grid_obj.y_total, grid_obj.z_total         = np.linspace(
            -0.5, +0.5, grid_obj.nx_total), np.linspace(-0.5, +0.5, grid_obj.ny_total), np.linspace(-0.5, +0.5, grid_obj.nz_total)
                                                                                                                    
        grid_points                                                  = np.vstack(np.meshgrid(grid_obj.x_total, grid_obj.y_total,
                                            grid_obj.z_total, indexing='ij')).reshape(3, -1).T
        cartesian_coord                                              = np.dot(grid_obj.A, grid_points.T).T
        grid_data                                                    = np.hstack((cartesian_coord, e_vals))
        f1                                                           = open(path_to_file, 'wb')
        np.savetxt(f4, grid_data)
        f1.close()

    # * Write the grid to vtk in unstructured points format
    def write_vts(grid_obj, path_to_file):
        """ Description
        :type grid_obj: instance of the grid3D class
        :param grid_obj: Contains all the info regarding the energy grid

        :type path_to_file: str
        :param path_to_file: path to the output file including full file name and extension.

        :raises:

        :rtype: Just writes the file.
        """

        import numpy as np
        from pyevtk.hl import pointsToVTK

        # For we need the grid x, y, z, E style
        # Define the crazy unstructured grid
        #
        dx, dy, dz                                                   = 1.0/grid_obj.nx_total, 1.0/grid_obj.ny_total, 1.0/grid_obj.nz_total
        X                                                            = np.arange(0, 1 + 0.1*dx, dx, dtype='float64')
        Y                                                            = np.arange(0, 1 + 0.1*dy, dy, dtype='float64')
        Z                                                            = np.arange(0, 1 + 0.1*dz, dz, dtype='float64')
        x                                                            = np.zeros((grid_obj.nx_total, grid_obj.ny_total, grid_obj.nz_total))
        y                                                            = np.zeros((grid_obj.nx_total, grid_obj.ny_total, grid_obj.nz_total))
        z                                                            = np.zeros((grid_obj.nx_total, grid_obj.ny_total, grid_obj.nz_total))
        for k in range(grid_obj.nz_total):
                for j in range(grid_obj.ny_total):
                        for i in range(grid_obj.nx_total):
                                x[i, j, k]                           = X[i]
                                y[i, j, k]                           = Y[j]
                                z[i, j, k]                           = Z[k]
        for k in range(grid_obj.nz_total):
                for j in range(grid_obj.ny_total):
                        for i in range(grid_obj.nx_total):
                                [x[i, j, k], y[i, j, k], z[i, j, k]] = np.dot(
                                    grid_obj.A, [x[i, j, k], y[i, j, k], z[i, j, k]])
        #Write pot into .vts file
        pointsToVTK(path_to_file, x, y, z, pointData={"Potential": grid_obj.pot_repeat})

    # * Converts the RASPA grid file into a vtk file, excluding the derivatives.
    def RASPA2vts(input_file, output_file):
        
    
        """ Description
        :type input_file: str   
        :param input_file: path to the RASPA .grid file
    
        :type output_file: str
        :param output_file:  path to the output file including full file name and extension. Filename should end with '.vts'   
        
        :raises:
    
        :rtype: Just writes the file.
        """

        import numpy as np
        import os
        import pandas as pd
        from pyevtk.hl import pointsToVTK

        data                                                         = pd.read_csv(input_file, delim_whitespace=True, header=None,names=['x', 'y', 'z', 'E', 'dEx', 'dEy', 'dEz'])
        

        x                                                            = data.x
        y                                                            = data.y
        z                                                            = data.z
        energy                                                       = data.E

        pointsToVTK(output_file, x, y, z, data={"energy": energy})
