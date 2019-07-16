import numpy as np
import pyIsoP.grid3D as grid3D 
import pyIsoP.potentials as potentials 
import pyIsoP.forcefields as forcefields 
import pyIsoP.writer as writer 

####################################################################
# Calculate the grid
t1=grid3D.grid3D('ZIF-4_mod.cif',spacing=0.5)
f1=forcefields.forcefields(t1,sigma=3.95, epsilon=46)
t2= grid3D.grid3D.grid_calc(t1,"lj",f1)

# Save coordinates for visualizing later
writer.writer.write_vts(t2,'zif-4_grid')
writer.writer.write_frame(t2,'zif-4_repeat.pdb')
#######################################################################


######################################################################
# Predict adsorption isotherm
############################################
import machlearn
import histo
import predictor

# Calculate/Read the histogram
h=histo.histo()
h=histo.histo.grid2histo(t2,h)

# Train the machine learning model (oie time thing)
ml=machlearn.machlearn(restarts=2)
ml=machlearn.machlearn.GPR4n1(ml,'for_gpr_large.dat',0.1)

# Isotherm Prediction
n_pressures=20
P=np.linspace(1E-5,100,n_pressures)
T=298
Vf=0.66
lcd=5.1
pld=2.4

X_test           = np.array([[Vf, lcd, pld],]*n_pressures)
X_test           = np.hstack((np.reshape(np.log10(P*1E5),(n_pressures,1)), X_test))
n1=ml.gp.predict(X_test)

# Predict isotherm
g_L_CH2=predictor.predictors.predict_isotherm(T,P,Vf,h,n1,epsilon=46,MA=14)


##########################################################################
#Automated molecular siting
#################################
m,e_min=grid3D.grid3D.find_minima(t2)
writer.writer.write_minima(t2, 'zif_minima.xyz')


# Link molecules to minima easily
import siter
s=siter.siter(t2,'Movie_ZIF-4_mod_2.2.2_303.000000_200000.000000_component_butane_0.pdb',4)
s=siter.siter.molecule2minima(s)
writer.writer.write_siter(s,'zif_4_siter.xyz')










