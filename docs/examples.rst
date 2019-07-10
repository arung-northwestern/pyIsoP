.. _examples:

===============================================================
Examples
===============================================================

Here we provide a few simple examples on using the different functionalities offered by PyIsoP.
Please refer to the corresponding API reference to find more about the modules and functions and the possible options.

.. _grid:

Energy Grid Calculation
=======================
PyIsoP uses a vectorized_ grid calculator in the :ref:`grid3D`, in conjunction with interatomic potentials in
:ref:`potentials` and parameters in RASPA format in the :ref:`forcefields` module.  

.. code-block:: python
       
        from grid3D import grid3D
        from potentials import potentials
        from forcefields import forcefields
        from writer import writer

        ####################################################################
        # Calculate the grid
        t1=grid3D('ZIF-4_mod.cif',spacing=0.5)          # Intialize grid3D object
        f1=forcefields(t1,sigma=3.95, epsilon=46)      # LJ parameters for the single-site probe 
        t2= grid3D.grid_calc(t1,"lj",f1)                          # Save to different object or overwrite the existing object

        # Save coordinates for visualizing later
        writer.writer.write_vts(t2,'zif-4_grid')                   # Write a binary vtk file
        writer.writer.write_frame(t2,'zif-4_repeat.pdb')    # Save the corresponding replicated structure corresponding to a 12.8 A (default) cut-off.

.. _pores:

Pore Structure Visualization
============================

Visualize the pores using the vtk file with any tool with volume rendering features to understand the pore
structure better. The image below is generated using Visit_ visualizer.

.. image:: ./images/zif-4.png
   :align: center
   :height: 600
   :width: 600

   

.. _histogram:

Energy Histogram
=========================================

PyIsoP contains the :ref:`histogram` module which offers 3 ways for the user to obtain the energy histogram.  The number of bins and :math:`E_max` can be set while initializing the histogram.
All the energies should be in the units of [K] to ensure consistency with the RASPA grid output.

1. From the PyIsoP :ref:`grid3D` object 

.. code-block:: python

        from histo import histo                 # import the histogram module
        h = histo()                                     # initialize a histo object
        h = histo.grid2histo(t2, h)            # update (overwrite) the histo object with histogram calculated from the grid3D object t2  

2.  Read in the energy grid from a RASPA style .grid file, with x, y, z, E data or from  .cube file. 

.. code-block:: python

        from histo import histo                 # import the histogram module
        h = histo()                                     # initialize a histo object
        h = histo.raspa2histo('raspa_grid_filename.grid' , ,h)            # update (overwrite) the histo object with histogram calculated from the RASPA grid file.
        h = histo.cube2histo('cube_filename.cube',h)            # update (overwrite) the histo object with histogram calculated from a .cube file

3. Read in the histogram as two column text file with no header. Bin-centers in one column,  normalized histogram in the other column.

.. code-block:: python

        from histo import histo                 # import the histogram module
        h = histo()                                     # initialize a histo object
        h = histo..file2histo('text_filename.dat', h)            # update (overwrite) the histo object with histogram calculated from the RASPA grid file.


.. _machlearn:

Coordination Number from Machine Learning
=========================================

In order to predict the guest-guest energy of hydrogen, we use a machine learning model (GPR) trained on the first-shell coordination number.
Please refer to :ref:`theory` section or our recent work by Gopalan *et al.*, :cite:`gopalan2019fast`  for more details. PyIsoP provides 
a pre-trained model at 77 K which can predict the hydrogen coordination numbers as a function of  [log10(P), void fraction, largest cavity diameter (A), pore limiting diameter (A)]

-   To load that model (details are in the  SI of the publication :cite:`gopalan2019fast` 

.. code-block:: python

        import joblib
        gp=joblib.load('gprmodel.joblib')               # Load the trained model
        n1 = gp.predict([logP, VF, LCD,PLD])           # Predict at 77 K for a set of  feature values for log10(pressure), void fraction, LCD and PLD in angstroms.

-   To train a new model using your own data (at your temperature of choice)  but with the default settings using Gaussian Process Regression, create a comma-separated-values (.csv) with 5 columns of "log(P)", "Vf", "lcd", "pld", "n1" with no header lines. Let's call it 'file_with_data.csv'

.. code-block:: python

        from machlearn import machlearn

        m1= machlearn(restarts=2)                          # Initialize object with  2 optimizer restarts
        m1 = machlearn.GPR4n1( m1, 'file_with_data.csv', 0.9)   # Train the model with 90 % training and 10 % Testing
        n1 = m1.predict([logP, VF, LCD,PLD])           # Predict at your temperature for a set of  feature values for log10(pressure), void fraction, LCD and PLD in angstroms.


-   Preferred:  To use algorithms other than GPR, users are encouraged to train their own model and be ready to provide :math:`n_1` as a vector (array corresponding to different pressures) to be fed into
    the :ref:`isotherm` calculation (example below) using the :ref:`predictor` .

.. _isotherm:

Adsorption Isotherm
==============================
PyIsoP takes in the temperature, pressures, void fraction, the energy histogram object, coordination numbers vector, Lennard-Jones well depth in [K] (should be consistent with the one used in the grid calculation) and the molecular weight (:math:`M_A`)
and predicts the adsorption isotherm in the units of grams per liter of the adsorbent. Combining all the examples from above, the isotherm can be calculated using the :ref:`predictor` as 

.. code-block:: python

        from predictor import predictors
        g_L_CH2=predictors.predict_isotherm(T,P,Vf,h,n1,epsilon=46,MA=14)


.. _screening:

High-throughput Screening 
===================
CoRE-MOF 2019 All Solvent Removed (12,914 structures)
-------------------------------------------------------------------

Using the same algorithm implemented as PyIsoP, we calculated the hydrogen adsorption isotherms for a
preliminary version of the CoRE MOF 2019-ASR (12,914 structures) from 1 Pa to 100 bar at two
temperatures (77 K and 160 K) in less than 24 hrs on 500 processors with a grid spacing of 1 :math:`\mathring{A}`.
The evolution of the gas uptake for the entire set of 12,914 materials at 77 K with increasing adsorption pressures
is depicted in the figures below. Having the entire isotherm enables us to answer
important questions regarding maximization of gas uptake quickly and accurately, like determining
the choice of the adsorption and desorption conditions for a material with given void fraction and
LCD or against any other textural property. For example, consider two materials A (highlighted in
blue in figures (e) and (f)) and B (highlighted in red in figures (e)
and (f)) with very similar void fractions, close to 0.85 but with different largest cavity
diameters of 13.5 :math:`\mathring{A}`  and 34.9 :math:`\mathring{A}`,  respectively. If one were to use A for storing hydrogen at 77
K, increasing the adsorption pressure from 42 bar \ref{fig:L5} to 100 bar \ref{fig:L6}) would give
an improvement of less than 1\% (56.939 g/L to 57.92 g/L) in gas uptake, hence is not worthwhile considering
the increased costs and risks of storing at higher pressures. Instead, if one were using B, the
same change in pressure will improve the the uptake by about 50\% (30.87 g/L to 44.065 g/L), which might be more
economically feasible. Please refer to Gopalan *et al.* :ref:`gopalan2019fast` for more information.



.. figure:: ./images/screening.png
  :align: center

.. _vectorized: https://numba.pydata.org/numba-doc/dev/user/vectorize.html
.. _VisIt: https://wci.llnl.gov/simulation/computer-codes/visit/