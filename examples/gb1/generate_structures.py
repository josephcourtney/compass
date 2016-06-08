#!/usr/bin/env python
import compass
import compass.mod
import compass.ros
import shutil, os
import modeller.automodel
protein='GB1'
compass.set_shiftx2_path( '/path/to/shiftx2-v107-mac-20120130/shiftx2.py' )
compass.set_working_directory( './working')
compass.set_model_directory('./models')
compass.set_sequence( 'MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE' )

if not os.path.exists(compass.model_directory):
    os.mkdir(compass.model_directory)
if not os.path.exists(compass.working_directory):
    os.mkdir(compass.working_directory)
shutil.copyfile(__file__, compass.model_directory+'/'+os.path.basename(__file__))

#choose one of the following sections

#To load generated models
# compass.set_pdb_regex('%s_\w{4}\.\w{9}\.pdb\Z' %protein)
# models = compass.load_pdbs()
# model_list=[model.pdb_path for model in models]

#To generate models from sequence
compass.mod.genali( compass.sequence, protein, compass.working_directory )
model_list = compass.mod.runmodeller(
	compass.working_directory + '/%s.ali' %protein,
	working=compass.working_directory,
        database_path = '/path/to/modeller_database',
	models_path = compass.model_directory,
	mod_per_temp=20,
	max_seq_id=0.98,
	min_seq_id = 0.15,
	num_iter = 6,
	max_eval = 0.8,
	gaps = True,
	loop_refine = False,
	md_level = modeller.automodel.refine.slow,
	repeat_optimization = 1
)


compass.ros.relaxlist( model_list, relaxprogram='/path/to/rosetta3.4/rosetta_source/bin/relax.default.macosgccrelease', rosettaDB='/path/to/rosetta3.4/rosetta_database', fast = False, CPUnum = 8 )
