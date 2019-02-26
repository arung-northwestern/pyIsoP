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
    TODO: Add function to write the vtk file
    TODO: Add function to write the energy values into a file.
    ! Check citation for xyz.py file 

        """
    import os
    import imp
    import numpy as np
    xyz_mod       = imp.load_source('xyz', 'xyz.py')  # A single file to make xyz file writing easy

    def __init__(self, *args, **kwargs):
        pass

    def write_details(grid_obj, path_to_folder):
        """ Description
        Reads in the grid3D object after calculation and writes out the important details into a named file
        :type grid_obj: An instance of the class grid3D
        :param grid_obj: grid3D object after calculation

        :type path_to_folder: str
        :param path_to_folder: Path to the folder where you want the file to be written to.

        :raises:

        :rtype:
        """
        stripped  = os.path.splitext(grid_obj.file)[0]
        f1        = open(stripped+'_details.dat', 'w')
        f1.write(str(grid_obj.nx_total)+'\t'+str(grid_obj.ny_total)+'\t'+str(grid_obj.nz_total)+'\n')
        f1.write(str(grid_obj.nx_cells)+'\t'+str(grid_obj.ny_cells)+'\t'+str(grid_obj.nz_cells)+'\n')
        f1.write(str(grid_obj.la)+'\t'+str(grid_obj.lb)+'\t'+str(grid_obj.lc)+'\n')
        f1.write(str(grid_obj.alpha)+'\t'+str(grid_obj.beta)+'\t'+str(grid_obj.gamma)+'\n')
        f1.close()

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
        stripped  = os.path.splitext(grid_obj.file)[0]
        f1        = open(stripped+'.xyz', 'w')
        out_coord = np.zeros((np.shape(grid_obj.coord)))
        out_coord = np.dot(A, grid_obj.coord.T).T
        xyz_mod.write_xyz(f1, out_coord, title=stripped+'.xyz', atomtypes=grid_obj.mof_atm_names)
        f1.close()


# def write_energy_values(grid_obj, path_to_file):
#     stripped    = os.path.splitext(grid_obj.file)[0]
#     f1          = open(stripped+'_E_values.dat', 'w')


# def write_vts(grid_obj, path_to_file):
