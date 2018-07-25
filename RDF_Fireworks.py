# This magnificent script calculates the RDF and the coordination numbers
# as fireworks, which can be kept track of easily


# Im still making up my mind on how to split the jobs
#  input creation
#  Running GCMC
#  Doing whatever after


from fireworks.core.firework import Firework, Workflow, Tracker
from fireworks.core.fworker import FWorker
from fireworks.core.launchpad import LaunchPad
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.user_objects.firetasks.script_task import ScriptTask
from fireworks.user_objects.firetasks.templatewriter_task import TemplateWriterTask
from fireworks.queue.queue_launcher import launch_rocket_to_queue #  rapidfire    
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter
import os
import pandas as pd
from pymatgen.io.cif import CifParser
import shutil


screening_folder   = '/home/agp971/SCOTTY/fireworks/screening'
cifs_folder        = '/projects/p20832/GCMC_DATA/Scotty_Metric/Hydrogen/1000coremofs_H2_1'
path_to_launchpad  = '/home/agp971/SCOTTY/fireworks/my_launchpad.yaml'
path_porediameters = screening_folder + '/coremof_small_set_pore_diameters.dat'
path_to_RDF_script = '/home/agp971/SCOTTY/RDF_quest.py'

# # setup the launchpad
# launchpad = LaunchPad(host='ds111258.mlab.com', port=11258, name='db_for_fw', username='arung', password='sivaliks')
launchpad = LaunchPad.from_file(path_to_launchpad)
#launchpad.reset('', require_password=False)


###################################################################################
# TASK 0: READ the specifications of the materials you want to set/query later
###################################################################################
PoreDiameters = pd.read_csv(path_porediameters, delim_whitespace=True)


##########################################################
# Calculate RDF as a fireworks
#########################################################


os.chdir(cifs_folder)
path  = os.getcwd()
folds = os.listdir('.')

for i in range(len(folds)):
    os.chdir(folds[i])
    #os.system('rm *.cif')
    folds1 = os.listdir('.')
    path2  = os.getcwd()

    for j in range(len(folds1)):
        os.chdir(folds1[j])
        path3       = os.getcwd()
        shutil.copy2(path_to_RDF_script, path3)
        task_string = 'cd ' + path3 + ' && python RDF_quest.py'
        # task1     = ScriptTask.from_str('export RASPA_DIR=/home/agp971/Pi_Complexation/RASPA-2.0/simulations')  # RASPA path
        task1       = ScriptTask.from_str(task_string)
        # task2     = ScriptTask.from_str('export RASPA_DIR=/home/agp971/Pi_Complexation/RASPA-2.0/simulations; /home/agp971/Pi_Complexation/RASPA-2.0/src/simulate $1')
        fw_name     = 'RDF_'+folds[i]+'_'+folds1[j]+'_Pa_'+'77_K'
        PLD         = PoreDiameters[PoreDiameters['Material'].str.contains(folds[i])]['PLD'].as_matrix()[0]
        LCD         = PoreDiameters[PoreDiameters['Material'].str.contains(folds[i])]['LCD'].as_matrix()[0]
        # tracker1  = Tracker(path3+'/Output/System_0/output_'+folds[i]+'.data')
        # tracker1  = Tracker(path3 + 'Output/System_0_0/' + os.listdir('Output/System_0')[0], nlines=20)
        submit_job  = Firework([task1], spec={"Pressure": float(folds1[j]), "Material": folds[i], "LCD": LCD, "PLD": PLD, "_trackers": []}, name=fw_name)  # Firetasks are always run sequentially.
        launchpad.add_wf(submit_job)
        os.chdir(path2)

    os.chdir(path)








# Test the rapidfire settings with the MOAB queue adapter.
# Launching the workflow into queue using the MOAB queue adapter
# moab = CommonAdapter('MOAB',q_name='normal', \
#      template_file=None,default_q_commands= \
#      {'MOAB': {'submit_cmd': 'msub', 'status_cmd': 'showq -w user=agp971'}})

# qadapter = CommonAdapter.from_file('my_qadapter.yaml')
# firetask = ScriptTask.from_str('echo "test fireworks queue" > howza.txt')
# test     = Firework(firetask)
# launchpad.add_wf(test)
# launchpad.add_wf(test)
# launchpad.add_wf(test)
# launchpad.add_wf(test)
# launch_rocket_to_queue(launchpad, FWorker(), qadapter,  \
#     launcher_dir='.', reserve=False, strm_lvl='INFO', create_launcher_dir=True, fill_mode=False, fw_id=None)
# rapidfire(launchpad, FWorker(), qadapter,  \
#     launch_dir='.', nlaunches='infinite', njobs_queue=2, reserve=False, strm_lvl='INFO', fill_mode=False)
# rapidf ire(launchpad, FWorker(), qadapter, launch_dir='.', njobs_queue=2, reserve=False, strm_lvl='INFO', fill_mode=False)
