#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import compass
import compass.mod
import compass.ros
import sys

protein='GB1'
compass.set_shiftx2_path( '/path/to/shiftx2-v107-mac-20120130/shiftx2.py' )
compass.set_working_directory( './working')
compass.set_model_directory('./models')
compass.set_sequence( 'MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE' )


ref_id = '2lgi'
compass_id = os.path.basename(os.path.abspath('.'))
print compass_id

compass.set_pdb_regex('GB1_\w{4}\.\w{9}_relaxed\.pdb\Z' )
reference_pdb = ( os.path.abspath('./reference/%s.pdb' % ref_id ) )

ref_model = compass.Model( reference_pdb )
ref_model.load_cosy( os.path.abspath( './reference/73gb1cosy.pks' ) )

shiftx_ref_model = compass.Model( reference_pdb )
try:
    shiftx_ref_model.pdb.main_struct.load_shifts( reference_pdb+'.shiftx' )
except IOError:
    compass.calculate_shifts( [shiftx_ref_model] )
shiftx_ref_model.cosy_peaks = compass.simulate_cosy( shiftx_ref_model )
compass.calculate_mhd_between_models( shiftx_ref_model, ref_model )
shiftx_ref_model.rmsd = compass.calculate_rmsd_between_models( shiftx_ref_model, ref_model )
models = compass.load_pdbs()


for model in models:
    try:
        model.pdb.main_struct.load_shifts( model.pdb_path+'.shiftx' )
    except IOError:
        compass.calculate_shifts( [model] )

for model in models:
    model.cosy_peaks = compass.simulate_cosy( model )
    compass.calculate_mhd_between_models( model, ref_model )
    compass.calculate_rmsd_between_models( model, ref_model )

out_file = open( os.path.dirname(os.path.realpath(__file__)) + "/results.txt", 'w')

for model in models+[shiftx_ref_model]:
    out_str = ''
    out_str += str(model.pdb_path) if model.pdb_path else '.'
    out_str += '\t'
    out_str += str(model.mhd) if model.mhd else '0.0'
    out_str += '\t'
    out_str += str(model.rmsd) if model.rmsd else '0.0'
    out_file.write(out_str+'\n')

rmat = compass.rmsd_mat( models )

out_file.write( '\n'.join( ['mhd\t'+('\t'.join([str(eb) for eb in ea])) for ea in compass.func_average_pairwise_rmsd( models, rmat, 'mhd' ) ] ) + '\n' )
out_file.write( '\n'.join( ['mhd_ref\t'+('\t'.join([str(eb) for eb in ea])) for ea in compass.func_average_pairwise_rmsd_vs_ref( models, ref_model, 'mhd' ) ] ) + '\n' )

