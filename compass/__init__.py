# -*- coding: utf-8 -*-

import numpy as np
import sys, urllib, re, math, shutil, os, string, random
import subprocess
from time import time
import compass.util

threetoone = {"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V", "TRP":"W", "TYR":"Y"}
onetothree = {v:k for k,v in threetoone.iteritems()}

# Settings

pdb_regex = None
working_directory = None
model_directory = None
shiftx2_path = None
sequence = None
cosy_peaks = None

class Logger:
    def __init__(self, name, tee = False ):
        self.logf = open(name, 'w', 1)
        self.stdout = sys.stdout
        sys.stdout = self
        self.tee = tee
    def close(self):
        sys.stdout = self.stdout
        self.logf.close()
    def write(self, data):
        self.logf.write(data)
        if self.tee:
            self.stdout.write(data)
    def __getattr__(self, name, *args):
        def method(*args):
            return getattr(self.stdout, name)
        return method(args)

class Model:
    def __init__(self, pdb_path):
        self.pdb_path = pdb_path
        self.pdb = PDB(open(pdb_path,'r').read())
        self.shiftx_path = ""
        self.cosy_path = ""
        self.cosy_peaks = []
        self.reference_pdb_path = ""
        self.reference_cosy_path = ""
        self.mhd = None
        self.rmsd = None
    def load_cosy(self, cosy_file ):
        f = open(cosy_file, 'r')
        self.cosy_peaks = []
        for line in f.readlines():
            self.cosy_peaks.append( [ float(e) for e  in line.split(',')[0:2] ] )
        return self.cosy_peaks

class Atom:
    def __init__(self, data_tup):
        self.pos = np.array( [ float(data_tup[8]), float(data_tup[9]), float(data_tup[10]) ] )
        self.mass = 0.0
        self.type = data_tup[2].strip()
        self.element = data_tup[13].strip()
        self.chain_id = data_tup[5].strip()
        self.res_type = data_tup[4].strip()
        self.res_num = int(data_tup[6])
        self.occupancy = float(data_tup[11])
        self.bfactor = float(data_tup[12])
        self.num = int(data_tup[1])
        self.model = 0
        self.chemical_shift = 0.0
    def pdb_str(self):
        return "%6s%5s  %4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s  "  % ('ATOM  ', self.num, str(self.type).ljust(4), self.res_type, self.chain_id, self.res_num, self.pos[0], self.pos[1], self.pos[2], self.occupancy, self.bfactor, str(self.element).rjust(12))
    @classmethod
    def from_string(cls, pdb_str):
      fields = [ (0,6), (6,11), (12,16), (16,17), (17,20), (21,22), (22,26), (26,27), (30,38), (38,46), (46,54), (54,60), (60,66), (72,76), (76,78), (78,80) ]
      return Atom( [ pdb_str[beg:end] for (beg,end) in fields ] )

class Residue:
    def __init__(self, res_num, res_type):
        self.index = res_num
        self.type = res_type
        self.atoms = dict()
    def add_atom(self, a):
        self.atoms[a.type] = a

class Structure:
    def __init__(self):
        self.residues = dict()
    def add_atom(self, a):
        if a.res_num not in self.residues.keys():
            self.residues[a.res_num] = Residue( a.res_num, a.res_type )
            self.residues[a.res_num].add_atom( a )
        else:
            self.residues[a.res_num].add_atom( a )
    def get_atom(self, res_num, type):
        return self.residues[res_num].atoms[type]
    def load_shifts(self, shiftx_file):
        shift_data = [ line.split() for line in open(shiftx_file,'r').readlines() ]
        for datum in shift_data:
            if int(datum[1]) in self.residues:
                res = self.residues[int(datum[1])]
                if datum[3] in res.atoms:
                    res.atoms[datum[3]].chemical_shift = float(datum[4])
    def pdb_str(self):
        s = ""
        for res in self.residues:
            for atom in res.atoms:
                s += atom.pdb_str()
        return s

class PDB:
    def __init__(self, pdb_str):
        current_chain = 0
        current_model = 0
        self.structures = dict()
        for line in pdb_str.split('\n'):
            rec_name = line[0:6]
            if rec_name == 'MODEL ':
                current_model = int(line.split()[1])
                self.structures[str(current_model)+str(current_chain)] = Structure()
            elif rec_name == 'ATOM  ':
                a = Atom.from_string( line )
                a.model = current_model
                if a.chain_id != current_chain:
                    current_chain = a.chain_id
                    self.structures[str(current_model)+a.chain_id] = Structure()
                self.structures[str(current_model)+a.chain_id].add_atom( a )
        struct_keys = [ (struct, len(self.structures[struct].residues) ) for struct in self.structures.keys() ]
        struct_keys.sort( key=lambda tup: tup[1] )
        self.main_struct = self.structures[struct_keys[-1][0]]

class Job:
    '''
        A Job object will be created for each subprocess
        models is a list of Models to be calculated for a subprocess
        popen is a subprocess.Popen object used to control a subprocess
        folder is the temp folder where you will store the pdb files to
    '''
    def __init__(self):
        self.models=[]
        self.popen=None
        self.folder=""

def set_sequence(str):
  global sequence
  sequence = str

def set_working_directory(str):
  global working_directory
  working_directory = str

def set_model_directory(str):
  global model_directory
  model_directory = str

def set_shiftx2_path(str):
  global shiftx2_path
  shiftx2_path = str

def set_pdb_regex(str):
  global pdb_regex
  pdb_regex = re.compile(str)

def load_peak_list(path, exp='COSY'):
  global cosy_peaks
  f = open(path, 'r')
  cosy_peaks = []
  for line in f.readlines():
    cosy_peaks.append( [ float(e) for e  in line.split(',')[0:2] ] )
  return cosy_peaks

def pdb_list():
  return [os.path.abspath(model_directory+'/'+f) for f in os.listdir(model_directory) if re.search(pdb_regex, f) ]

def load_pdbs( pdb_files = None ):
    '''Returns a list of TestModel objects given a list of PDB files'''
    if pdb_files is None:
      return [ Model( f ) for f in pdb_list() ]
    else:
      return [ Model( f ) for f in pdb_files]

def calculate_shifts( model_list, shiftx2=None, CPUnum=8, working='./', replace = True, logger = False ):
    '''acts as an interfase for shiftx.
       Input: a list of pdb files
       Output: a list of shiftx output files
       Uses subprocess.Popen to create subprocesses. Shiftx2 is run in batch mode. Temp file folders will
       be created in the working folder. Any two models cannot share the same pdb filename.
    '''

    print "Calculating chemical shifts with SHIFTX2 parallelized into %d threads for:" % CPUnum
    print '\n'.join( [ model.pdb_path for model in model_list ] ), '\n'

    if type(model_list) is not list:
      model_list = [model_list]

    if shiftx2 is None:
      if shiftx2_path is None:
        shiftx2=os.environ['COMPASSInstall']+'/shiftx2-v107-mac-20120130/shiftx2.py'
      else:
        shiftx2 = shiftx2_path

    shiftx_files = []

    #divide jobs
    #processlist is a list of Job objects
    processlist=[]
    for i in range(CPUnum):
        processlist.append(Job())
    counter=0
    for model in model_list:
        if counter==CPUnum:
            counter=0
        model.shiftx_path = ( "%s.shiftx" % model.pdb_path )
        shiftx_files.append( model.shiftx_path )
        processlist[counter].models.append(model)
        counter+=1

    #in case you have less than CPUnum of pdb files
    for i in range(len(processlist)):
        if processlist[i].models==[]:
            processlist=processlist[:i]
            break

    #create a unique list of temporary directory names for the form temp********** where * is a random alphanumeric character
    characters = string.lowercase+string.uppercase+string.digits
    temp_dir_names = []
    while len(temp_dir_names) < len(processlist):
      new_temp_dir = working+'temp'+''.join(random.sample(characters,10))
      if new_temp_dir not in temp_dir_names:
        temp_dir_names.append( new_temp_dir )

    #create the temp folders and copy the pdb files
    for i in range(len(processlist)):
        processlist[i].folder = temp_dir_names[i]
        if os.path.exists(processlist[i].folder):
            shutil.rmtree(processlist[i].folder)
        os.mkdir(processlist[i].folder)
        processlist[i].folder=os.path.abspath(processlist[i].folder)
        for model in processlist[i].models:
            shutil.copy(model.pdb_path,processlist[i].folder)
    initial_path=os.getcwd()
    #create calculation subprocesses
    for i in range(len(processlist)):
        os.chdir(processlist[i].folder)
        argument = ['python', shiftx2, '-b', '*.pdb', '-f', 'BMRB', '-a', 'ALL']
        if logger:
            processlist[i].popen = subprocess.Popen( args = argument, stderr = subprocess.STDOUT, stdout = subprocess.PIPE )
        else:
            processlist[i].popen = subprocess.Popen( args = argument )
        print ('Command for process %d:\t' % processlist[i].popen.pid) + ' '.join(argument)

    os.chdir(initial_path)
    #wait calculations to finish
    for process in processlist:
        if process.popen.poll()==None:
            process.popen.wait()
    #rename and move the '*.cs' files to model.shiftx_path
    #remove the temp folders
    for process in processlist:
        for model in process.models:
            scr = process.folder+'/'+os.path.basename(model.pdb_path)+'.cs'
            os.rename(scr,model.pdb_path+".shiftx")
        shutil.rmtree(process.folder)
        if logger:
            print "####### output for process %d #######" % process.popen.pid
            print process.popen.stdout.read()


    if replace:
      for model in model_list:
        model.pdb.main_struct.load_shifts(model.shiftx_path)

    return shiftx_files

def simulate_cosy( model ):
  '''simulates a cosy file for a list of chemical shifts
     Input: a shiftx output file path
     Output: a filename for a cosy peak list file
  '''
  aa_atom_names = {
    'ALA': ['N','CA','CB','C'],
    'CYS': ['N','CA','CB','C'],
    'ASP': ['N','CA','CB','CG','C'],
    'GLU': ['N','CA','CB','CG','CD','C'],
    'PHE': ['N','CA','CB','CG','CD1','CD2','CE1','CE2','CZ','C'],
    'GLY': ['N','CA','C'],
    'HIS': ['N','CA','CB','CG','ND1','CD2','CE1','NE2','C'],
    'ILE': ['N','CA','CB','CG1','CG2','CD1','C'],
    'LYS': ['N','CA','CB','CG','CD','CE','NZ','C'],
    'LEU': ['N','CA','CB','CG','CD1','CD2','C'],
    'MET': ['N','CA','CB','CG','CE','C'],
    'ASN': ['N','CA','CB','CG','ND','C'],
    'PRO': ['N','CA','CB','CG','CD','C'],
    'GLN': ['N','CA','CB','CG','CD','NE2', 'C'],
    'ARG': ['N','CA','CB','CG','CD','NE','CZ','NH1','NH2','C'],
    'SER': ['N','CA','CB','C'],
    'THR': ['N','CA','CB','CG2','C'],
    'VAL': ['N','CA','CB','CG1','CG2','C'],
    'TRP': ['N','CA','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2','C'],
    'TYR': ['N','CA','CB','CG','CD1','CD2','CE1','CE2','CZ','C']
  }

  aa_dist_graphs = {
    'ALA': [
        [0,1,2,2],
        [1,0,1,1],
        [2,1,0,2],
        [2,1,2,0]
    ],
    'CYS': [
        [0,1,2,2],
        [1,0,1,1],
        [2,1,0,2],
        [2,1,2,0]
    ],
    'ASP': [
        [0,1,2,3,2],
        [1,0,1,2,1],
        [2,1,0,1,2],
        [3,2,1,0,3],
        [2,1,2,3,0]
    ],
    'GLU': [
        [0,1,2,3,4,2],
        [1,0,1,2,3,1],
        [2,1,0,1,2,2],
        [3,2,1,0,1,3],
        [4,3,2,1,0,4],
        [2,1,2,3,4,0]
    ],
    'PHE': [
        [0,1,2,3,4,4,5,5,6,2],
        [1,0,1,2,3,3,4,4,5,1],
        [2,1,0,1,2,2,3,3,4,2],
        [3,2,1,0,1,1,2,2,3,3],
        [4,3,2,1,0,2,1,3,2,4],
        [4,3,2,1,2,0,3,1,2,4],
        [5,4,3,2,1,3,0,2,1,5],
        [5,4,3,2,3,1,2,0,1,5],
        [6,5,4,3,2,2,1,1,0,6],
        [2,1,2,3,4,4,5,5,6,0]
    ],
    'GLY': [
        [0,1,2],
        [1,0,1],
        [2,1,0]
    ],
    'HIS': [
        [0,1,2,3,4,4,5,5,2],
        [1,0,1,2,3,3,4,4,1],
        [2,1,0,1,2,2,3,3,2],
        [3,2,1,0,1,1,2,2,3],
        [4,3,2,1,0,2,1,2,4],
        [4,3,2,1,2,0,2,1,4],
        [5,4,3,2,1,2,0,1,5],
        [5,4,3,2,2,1,1,0,5],
        [2,1,2,3,4,4,5,5,0]
    ],
    'ILE': [
        [0,1,2,3,3,4,2],
        [1,0,1,2,2,3,1],
        [2,1,0,1,1,2,2],
        [3,2,1,0,2,1,3],
        [3,2,1,2,0,3,3],
        [4,3,2,1,3,0,4],
        [2,1,2,3,3,4,0]
    ],
    'LYS': [
        [0,1,2,3,4,5,6,2],
        [1,0,1,2,3,4,5,1],
        [2,1,0,1,2,3,4,2],
        [3,2,1,0,1,2,3,3],
        [4,3,2,1,0,1,2,4],
        [5,4,3,2,1,0,1,5],
        [6,5,4,3,2,1,0,6],
        [2,1,2,3,4,5,6,0]
    ],
    'LEU': [
        [0,1,2,3,4,4,2],
        [1,0,1,2,3,3,1],
        [2,1,0,1,2,2,2],
        [3,2,1,0,1,1,3],
        [4,3,2,1,0,2,4],
        [4,3,2,1,2,0,4],
        [2,1,2,3,4,4,0]
    ],
    'MET': [
        [0,1,2,3,5,2],
        [1,0,1,2,4,1],
        [2,1,0,1,3,2],
        [3,2,1,0,2,3],
        [5,4,3,2,0,5],
        [2,1,2,3,5,0]
    ],
    'ASN': [
        [0,1,2,3,4,2],
        [1,0,1,2,3,1],
        [2,1,0,1,2,2],
        [3,2,1,0,1,3],
        [4,3,2,1,0,4],
        [2,1,2,3,4,0]
    ],
    'PRO': [
        [0,1,2,2,1,2],
        [1,0,1,2,2,1],
        [2,1,0,1,2,2],
        [2,2,1,0,1,3],
        [1,2,2,1,0,3],
        [2,1,2,3,3,0]
    ],
    'GLN': [
        [0,1,2,3,4,5,2],
        [1,0,1,2,3,4,1],
        [2,1,0,1,2,3,2],
        [3,2,1,0,1,2,3],
        [4,3,2,1,0,1,4],
        [5,4,3,2,1,0,5],
        [2,1,2,3,4,5,0]
    ],
    'ARG': [
        [0,1,2,3,4,5,6,7,7,2],
        [1,0,1,2,3,4,5,6,6,1],
        [2,1,0,1,2,3,4,5,5,2],
        [3,2,1,0,1,2,3,4,4,3],
        [4,3,2,1,0,1,2,3,3,4],
        [5,4,3,2,1,0,1,2,2,5],
        [6,5,4,3,2,1,0,1,1,6],
        [7,6,5,4,3,2,1,0,2,7],
        [7,6,5,4,3,2,1,2,0,7],
        [2,1,2,3,4,5,6,7,7,0]
    ],
    'SER': [
        [0,1,2,2],
        [1,0,1,1],
        [2,1,0,2],
        [2,1,2,0]
    ],
    'THR': [
        [0,1,2,3,2],
        [1,0,1,2,1],
        [2,1,0,1,2],
        [3,2,1,0,3],
        [2,1,2,3,0]
    ],
    'VAL': [
        [0,1,2,3,3,2],
        [1,0,1,2,2,1],
        [2,1,0,1,1,2],
        [3,2,1,0,2,3],
        [3,2,1,2,0,3],
        [2,1,2,3,3,0]
    ],
    'TRP': [
        [0,1,2,3,4,4,5,5,5,6,6,7,2],
        [1,0,1,2,3,3,4,4,4,5,5,6,1],
        [2,1,0,1,2,2,3,3,3,4,4,5,2],
        [3,2,1,0,1,1,2,2,2,3,3,4,3],
        [4,3,2,1,0,2,1,2,3,3,4,4,4],
        [4,3,2,1,2,0,2,1,1,2,2,3,4],
        [5,4,3,2,1,2,0,1,3,2,4,3,5],
        [5,4,3,2,2,1,1,0,2,1,3,2,5],
        [5,4,3,2,3,1,3,2,0,3,1,2,5],
        [6,5,4,3,3,2,2,1,3,0,2,1,6],
        [6,5,4,3,4,2,4,3,1,2,0,1,6],
        [7,6,5,4,4,3,3,2,2,1,1,0,7],
        [2,1,2,3,4,4,5,5,5,6,6,7,0]
    ],
    'TYR': [
        [0,1,2,3,4,4,5,5,6,2],
        [1,0,1,2,3,3,4,4,5,1],
        [2,1,0,1,2,2,3,3,4,2],
        [3,2,1,0,1,1,2,2,3,3],
        [4,3,2,1,0,2,1,3,2,4],
        [4,3,2,1,2,0,3,1,2,4],
        [5,4,3,2,1,3,0,2,1,5],
        [5,4,3,2,3,1,2,0,1,5],
        [6,5,4,3,2,2,1,1,0,6],
        [2,1,2,3,4,4,5,5,6,0]
    ],
  }

  cosy_peaks = []

  struct_keys = [ (struct, len(model.pdb.structures[struct].residues) ) for struct in model.pdb.structures.keys() ]
  struct_keys.sort( key=lambda tup: tup[1] )

  structure = model.pdb.structures[struct_keys[-1][0]]

  for res in structure.residues.keys():
    for atom_a in structure.residues[res].atoms.keys():
      res_type = structure.residues[res].type
      for atom_b in structure.residues[res].atoms.keys():
        if atom_a[0]=='C' and atom_b[0]=='C' and atom_a in aa_atom_names[res_type] and atom_b in aa_atom_names[res_type]:
          ia = aa_atom_names[res_type].index(atom_a)
          ib = aa_atom_names[res_type].index(atom_b)
          if aa_dist_graphs[res_type][ia][ib] == 1:
            cosy_peaks.append( (structure.residues[res].atoms[atom_a].chemical_shift, structure.residues[res].atoms[atom_b].chemical_shift ) )
  return cosy_peaks

def simulate_cosy_with_label( model ):
  '''simulates a cosy file for a list of chemical shifts
     Input: a shiftx output file path
     Output: a filename for a cosy peak list file
  '''
  aa_atom_names = {
    'ALA': ['N','CA','CB','C'],
    'CYS': ['N','CA','CB','C'],
    'ASP': ['N','CA','CB','CG','C'],
    'GLU': ['N','CA','CB','CG','CD','C'],
    'PHE': ['N','CA','CB','CG','CD1','CD2','CE1','CE2','CZ','C'],
    'GLY': ['N','CA','C'],
    'HIS': ['N','CA','CB','CG','ND1','CD2','CE1','NE2','C'],
    'ILE': ['N','CA','CB','CG1','CG2','CD1','C'],
    'LYS': ['N','CA','CB','CG','CD','CE','NZ','C'],
    'LEU': ['N','CA','CB','CG','CD1','CD2','C'],
    'MET': ['N','CA','CB','CG','CE','C'],
    'ASN': ['N','CA','CB','CG','ND','C'],
    'PRO': ['N','CA','CB','CG','CD','C'],
    'GLN': ['N','CA','CB','CG','CD','NE2', 'C'],
    'ARG': ['N','CA','CB','CG','CD','NE','CZ','NH1','NH2','C'],
    'SER': ['N','CA','CB','C'],
    'THR': ['N','CA','CB','CG2','C'],
    'VAL': ['N','CA','CB','CG1','CG2','C'],
    'TRP': ['N','CA','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2','C'],
    'TYR': ['N','CA','CB','CG','CD1','CD2','CE1','CE2','CZ','C']
  }

  aa_dist_graphs = {
    'ALA': [
        [0,1,2,2],
        [1,0,1,1],
        [2,1,0,2],
        [2,1,2,0]
    ],
    'CYS': [
        [0,1,2,2],
        [1,0,1,1],
        [2,1,0,2],
        [2,1,2,0]
    ],
    'ASP': [
        [0,1,2,3,2],
        [1,0,1,2,1],
        [2,1,0,1,2],
        [3,2,1,0,3],
        [2,1,2,3,0]
    ],
    'GLU': [
        [0,1,2,3,4,2],
        [1,0,1,2,3,1],
        [2,1,0,1,2,2],
        [3,2,1,0,1,3],
        [4,3,2,1,0,4],
        [2,1,2,3,4,0]
    ],
    'PHE': [
        [0,1,2,3,4,4,5,5,6,2],
        [1,0,1,2,3,3,4,4,5,1],
        [2,1,0,1,2,2,3,3,4,2],
        [3,2,1,0,1,1,2,2,3,3],
        [4,3,2,1,0,2,1,3,2,4],
        [4,3,2,1,2,0,3,1,2,4],
        [5,4,3,2,1,3,0,2,1,5],
        [5,4,3,2,3,1,2,0,1,5],
        [6,5,4,3,2,2,1,1,0,6],
        [2,1,2,3,4,4,5,5,6,0]
    ],
    'GLY': [
        [0,1,2],
        [1,0,1],
        [2,1,0]
    ],
    'HIS': [
        [0,1,2,3,4,4,5,5,2],
        [1,0,1,2,3,3,4,4,1],
        [2,1,0,1,2,2,3,3,2],
        [3,2,1,0,1,1,2,2,3],
        [4,3,2,1,0,2,1,2,4],
        [4,3,2,1,2,0,2,1,4],
        [5,4,3,2,1,2,0,1,5],
        [5,4,3,2,2,1,1,0,5],
        [2,1,2,3,4,4,5,5,0]
    ],
    'ILE': [
        [0,1,2,3,3,4,2],
        [1,0,1,2,2,3,1],
        [2,1,0,1,1,2,2],
        [3,2,1,0,2,1,3],
        [3,2,1,2,0,3,3],
        [4,3,2,1,3,0,4],
        [2,1,2,3,3,4,0]
    ],
    'LYS': [
        [0,1,2,3,4,5,6,2],
        [1,0,1,2,3,4,5,1],
        [2,1,0,1,2,3,4,2],
        [3,2,1,0,1,2,3,3],
        [4,3,2,1,0,1,2,4],
        [5,4,3,2,1,0,1,5],
        [6,5,4,3,2,1,0,6],
        [2,1,2,3,4,5,6,0]
    ],
    'LEU': [
        [0,1,2,3,4,4,2],
        [1,0,1,2,3,3,1],
        [2,1,0,1,2,2,2],
        [3,2,1,0,1,1,3],
        [4,3,2,1,0,2,4],
        [4,3,2,1,2,0,4],
        [2,1,2,3,4,4,0]
    ],
    'MET': [
        [0,1,2,3,5,2],
        [1,0,1,2,4,1],
        [2,1,0,1,3,2],
        [3,2,1,0,2,3],
        [5,4,3,2,0,5],
        [2,1,2,3,5,0]
    ],
    'ASN': [
        [0,1,2,3,4,2],
        [1,0,1,2,3,1],
        [2,1,0,1,2,2],
        [3,2,1,0,1,3],
        [4,3,2,1,0,4],
        [2,1,2,3,4,0]
    ],
    'PRO': [
        [0,1,2,2,1,2],
        [1,0,1,2,2,1],
        [2,1,0,1,2,2],
        [2,2,1,0,1,3],
        [1,2,2,1,0,3],
        [2,1,2,3,3,0]
    ],
    'GLN': [
        [0,1,2,3,4,5,2],
        [1,0,1,2,3,4,1],
        [2,1,0,1,2,3,2],
        [3,2,1,0,1,2,3],
        [4,3,2,1,0,1,4],
        [5,4,3,2,1,0,5],
        [2,1,2,3,4,5,0]
    ],
    'ARG': [
        [0,1,2,3,4,5,6,7,7,2],
        [1,0,1,2,3,4,5,6,6,1],
        [2,1,0,1,2,3,4,5,5,2],
        [3,2,1,0,1,2,3,4,4,3],
        [4,3,2,1,0,1,2,3,3,4],
        [5,4,3,2,1,0,1,2,2,5],
        [6,5,4,3,2,1,0,1,1,6],
        [7,6,5,4,3,2,1,0,2,7],
        [7,6,5,4,3,2,1,2,0,7],
        [2,1,2,3,4,5,6,7,7,0]
    ],
    'SER': [
        [0,1,2,2],
        [1,0,1,1],
        [2,1,0,2],
        [2,1,2,0]
    ],
    'THR': [
        [0,1,2,3,2],
        [1,0,1,2,1],
        [2,1,0,1,2],
        [3,2,1,0,3],
        [2,1,2,3,0]
    ],
    'VAL': [
        [0,1,2,3,3,2],
        [1,0,1,2,2,1],
        [2,1,0,1,1,2],
        [3,2,1,0,2,3],
        [3,2,1,2,0,3],
        [2,1,2,3,3,0]
    ],
    'TRP': [
        [0,1,2,3,4,4,5,5,5,6,6,7,2],
        [1,0,1,2,3,3,4,4,4,5,5,6,1],
        [2,1,0,1,2,2,3,3,3,4,4,5,2],
        [3,2,1,0,1,1,2,2,2,3,3,4,3],
        [4,3,2,1,0,2,1,2,3,3,4,4,4],
        [4,3,2,1,2,0,2,1,1,2,2,3,4],
        [5,4,3,2,1,2,0,1,3,2,4,3,5],
        [5,4,3,2,2,1,1,0,2,1,3,2,5],
        [5,4,3,2,3,1,3,2,0,3,1,2,5],
        [6,5,4,3,3,2,2,1,3,0,2,1,6],
        [6,5,4,3,4,2,4,3,1,2,0,1,6],
        [7,6,5,4,4,3,3,2,2,1,1,0,7],
        [2,1,2,3,4,4,5,5,5,6,6,7,0]
    ],
    'TYR': [
        [0,1,2,3,4,4,5,5,6,2],
        [1,0,1,2,3,3,4,4,5,1],
        [2,1,0,1,2,2,3,3,4,2],
        [3,2,1,0,1,1,2,2,3,3],
        [4,3,2,1,0,2,1,3,2,4],
        [4,3,2,1,2,0,3,1,2,4],
        [5,4,3,2,1,3,0,2,1,5],
        [5,4,3,2,3,1,2,0,1,5],
        [6,5,4,3,2,2,1,1,0,6],
        [2,1,2,3,4,4,5,5,6,0]
    ],
  }

  cosy_peaks = []

  struct_keys = [ (struct, len(model.pdb.structures[struct].residues) ) for struct in model.pdb.structures.keys() ]
  struct_keys.sort( key=lambda tup: tup[1] )

  structure = model.pdb.structures[struct_keys[-1][0]]

  for res in structure.residues.keys():
    for atom_a in structure.residues[res].atoms.keys():
      res_type = structure.residues[res].type
      for atom_b in structure.residues[res].atoms.keys():
        if atom_a[0]=='C' and atom_b[0]=='C' and atom_a in aa_atom_names[res_type] and atom_b in aa_atom_names[res_type]:
          ia = aa_atom_names[res_type].index(atom_a)
          ib = aa_atom_names[res_type].index(atom_b)
          if aa_dist_graphs[res_type][ia][ib] == 1:
            cosy_peaks.append( ((structure.residues[res].atoms[atom_a].chemical_shift, structure.residues[res].atoms[atom_b].chemical_shift), (structure.residues[res].atoms[atom_a].res_type, structure.residues[res].atoms[atom_a].type), (structure.residues[res].atoms[atom_b].res_type, structure.residues[res].atoms[atom_b].type) ) )
  return cosy_peaks


def calculate_mhd( model, cs_min = 10, cs_max = 80, diagonal_tol = 0.0 ):
    '''calculates the modified Hausdorff distance between two cosy spectra
       Input: a test cosy peak list file and a reference cosy peak list file
       Output: the modified Hausdorff distance as a float
    '''
    if not cosy_peaks or len(cosy_peaks) < 1:
      raise "Experimental COSY peak list is empty! Set it with compass.load_peak_list()"
    model.reference_pdb_path = None

    pruned_model_peaks = [ peak for peak in model.cosy_peaks if min(peak) > cs_min and max(peak) < cs_max and ( math.abs(peak[0]-peak[1]) / 1.414213562373095 ) > diagonal_tol ]
    pruned_exp_peaks = [ peak for peak in cosy_peaks if min(peak) > cs_min and max(peak) < cs_max and ( math.abs(peak[0]-peak[1]) / 1.414213562373095 ) > diagonal_tol ]

    model.mhd = compass.util.mhd( pruned_model_peaks, pruned_exp_peaks )
    return model.mhd

def calculate_mhd_between_models( model, reference_model, cs_min = 10, cs_max = 80, diagonal_tol = 0.0 ):
    '''calculates the modified Hausdorff distance between two cosy spectra
       Input: a test cosy peak list file and a reference cosy peak list file
       Output: the modified Hausdorff distance as a float
    '''
    model.reference_pdb_path = reference_model.pdb_path
    model.reference_cosy_path = reference_model.cosy_path

    pruned_model_peaks = [ peak for peak in model.cosy_peaks if min(peak) > cs_min and max(peak) < cs_max and ( abs(peak[0]-peak[1]) / 1.414213562373095 ) > diagonal_tol ]
    pruned_ref_peaks = [ peak for peak in reference_model.cosy_peaks if min(peak) > cs_min and max(peak) < cs_max and ( abs(peak[0]-peak[1]) / 1.414213562373095 ) > diagonal_tol ]

    model.mhd = compass.util.mhd( pruned_model_peaks, pruned_ref_peaks )
    return model.mhd

def calculate_rmsd_between_models( model, reference_model ):
    '''calculates the root mean square distance between two pdb files
       Input: a test pdb file and a reference pdb file
       Output: the rmsd as a float
    '''
    model.reference_pdb_path = reference_model.pdb_path
    model.reference_cosy_path = reference_model.cosy_path

    model.rmsd = ca_rmsd( model, reference_model )
    return model.rmsd

def generate_report( models, ref_model, name ):
    '''generates a report for the calculation'''
    buf = ""
    buf += "This is a report"
    f = open(name, 'w')
    f.write(buf)
    return buf

# RMSD code

def fetch_pdb_as_string(id):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    return urllib.urlopen(url).read()

def sq_dist( pos_1, pos_2 ):
    return (pos_1[0] - pos_2[0])**2 + (pos_1[1] - pos_2[1])**2 + (pos_1[2] - pos_2[2])**2

def ca_rmsd( model_a, model_b ):
    '''Calcualted the RMSD between the alpha carbons or model 0 chain A of two PDB objects using the Kabsh algorith'''

    def get_seq(model):
        str = ''
        for i in range(min(model.pdb.main_struct.residues.keys()),(max(model.pdb.main_struct.residues.keys())+1)):
            if i in model.pdb.main_struct.residues:
                str += threetoone[model.pdb.main_struct.residues[i].type]
            else:
                str += '-'
        return str

    def get_cas(model):
        cas = []
        for i in range(min(model.pdb.main_struct.residues.keys()),(max(model.pdb.main_struct.residues.keys())+1)):
            if i in model.pdb.main_struct.residues:
                if 'CA' in model.pdb.main_struct.residues[i].atoms:
                    cas.append( model.pdb.main_struct.residues[i].atoms['CA'] )
                else:
                    cas.append( None )
            else:
                cas.append( None )
        return cas

    def find_shift( seq_a, seq_b ):
        mat = np.zeros( (len(seq_a), len(seq_b)) )
        minlen = min(len(seq_a),len(seq_b))
        for ia, ea in enumerate(seq_a):
            for ib, eb in enumerate(seq_b):
                if ea == eb and ea != '-':
                    mat[ia,ib] = 0 + abs(float(ia-ib)/minlen)
                else:
                    mat[ia,ib] = 1 + abs(float(ia-ib)/minlen)
        trace_sum = [ (np.trace(mat, offset=i) / len(np.diag(mat, k=i)),i) for i in range(-len(seq_a)+1,len(seq_b)) ]
        return -min(trace_sum)[1]

    def aligned_subseqs( model_a, seq_b ):
        shift = find_shift(get_seq(model_a), get_seq(model_b))
        filt_seq_a = get_cas(model_a)[shift:]
        filt_seq_b = get_cas(model_b)[0:len(filt_seq_a)]
        l = min(len(filt_seq_a), len(filt_seq_b))
        return filt_seq_a[0:l], filt_seq_b[0:l]

    pa_cas, pb_cas = aligned_subseqs( model_a, model_b )

    pa_num_cas = len(list(pa_cas))
    pb_num_cas = len(list(pb_cas))

    pa_centroid = reduce( lambda a, b: a+b, (atom.pos for atom in pa_cas) ) / pa_num_cas
    pb_centroid = reduce( lambda a, b: a+b, (atom.pos for atom in pb_cas) ) / pb_num_cas

    pa_pos_mat = np.array( [ a.pos - pa_centroid for a in pa_cas ] )
    pb_pos_mat = np.array( [ a.pos - pb_centroid for a in pb_cas ] )

    cov_mat = np.dot( pa_pos_mat.transpose(), pb_pos_mat )
    U, s, Vt = np.linalg.svd( cov_mat )

    r = np.identity(3)
    r[2,2] = np.sign( np.linalg.det( np.dot( Vt.transpose(), U.transpose() ) ) )
    s[-1] = s[-1] * r[2,2]

    R = np.dot( np.dot( Vt.transpose(), r ), U.transpose() )

    msd = ( ( sum(sum(pa_pos_mat * pa_pos_mat)) + sum(sum(pb_pos_mat * pb_pos_mat)) - 2*sum(s) ) / pa_num_cas )
    if msd < 0:
        rmsd = 0.0
    else:
        rmsd = np.sqrt(msd)

    return rmsd

def average_pairwise_rmsd( model ):
    sum = 0.0
    tot_pairs = (len(models)*(len(models)+1)) / 2 - 1
    
    for ia,ea in enumerate(models):
        for eb in models[0,ia]:
            sum += compass.calculate_rmsd_between_models(ea, eb)
    return sum / tot_pairs

def rmsd_mat(models):
    rmsd_mat = np.zeros((len(models), len(models)))
    for ia,ea in enumerate(models):
        for ib,eb in enumerate(models):
            if rmsd_mat[ia,ib] == 0:
                rmsd_mat[ia,ib] = compass.ca_rmsd(ea,eb)
    return rmsd_mat

def func_average_pairwise_rmsd(models, rmat, var):
    assert len(models) != 0
    ordered_indices = zip(*sorted( list(enumerate(models)), key=lambda m: getattr(m[1], var) ))[0]
    sum = 0.0
    num = 0.0
    tot_pairs = (len(models)*(len(models)+1)) / 2 - 1
    
    rmsd_list = []
    for ia in ordered_indices:
        for ib in ordered_indices[0:ia]:
            sum += rmat[ia,ib]
            num += 1
        if num > 0:
            rmsd_list.append(((num/tot_pairs), (sum/num)))
    return rmsd_list

def func_average_pairwise_rmsd_vs_ref(models, ref_model, var):
    assert len(models) != 0
    ordered_indices = zip(*sorted( list(enumerate(models)), key=lambda m: getattr(m[1], var) ))[0]
    sum = 0.0
    num = 0.0
    tot_pairs = (len(models)) / 2 - 1
    
    rmsd_list = []
    for ia in ordered_indices:
        sum += compass.calculate_rmsd_between_models(models[ia], ref_model)
        num += 1
        if num > 0:
            rmsd_list.append(((num/tot_pairs), (sum/num)))
    return rmsd_list

def strtoshiftx(input_file,output_file='default',format='default'):
    '''
    A function to convert .str file to .shiftx file
    Inputs:
        input_file: A path to the .str file
        output_file: A path to the .shiftx file. When set to 'default', the shiftx file will be saved with the same name of the input file with the extension .shiftx
        format: format string indicating the format of the lines in the .shfitx file.
    Function return:
        A path to the output file   
    '''
    infilename=os.path.abspath(input_file)
    if output_file=='default':
        if infilename.endswith('.str'):
            outfilename=infilename[:-4]+'.shiftx'
        else:
            outfilename=infilename+'.shiftx'
    else:
        outfilename=os.path.abspath(output_file)

    infile=open(infilename,'r')
    outfile=open(outfilename,'w')
    
    columns=[]
    data=[]
    for line in infile:
        linelist=line.split()
        if linelist==[]:
            continue
        if columns==[] and not'_Atom_chem_shift.' in linelist[0]:
            continue
        if 'stop_' in linelist and data!=[]:
            break
        if '_Atom_chem_shift.' in linelist[0]:
            columns.append(linelist[0])
            continue
        if len(linelist)==len(columns) and len(columns)>1:
            data.append(linelist)
    IDindex=columns.index('_Atom_chem_shift.ID')
    SeqIDindex=columns.index('_Atom_chem_shift.Seq_ID')
    Comp_IDindex=columns.index('_Atom_chem_shift.Comp_ID')
    Atom_IDindex=columns.index('_Atom_chem_shift.Atom_ID')
    Valindex=columns.index('_Atom_chem_shift.Val')
    IDlen=str(len(data[-1][IDindex])+1)
    SeqIDlen=str(len(data[-1][SeqIDindex])+1)
    Comp_IDlen=str(max([len(line[Comp_IDindex]) for line in data])+1)
    Atom_IDlen=str(max([len(line[Atom_IDindex]) for line in data])+1)
    Vallen=str(max([len(line[Valindex]) for line in data])+1)
    
    if format=='default':
        formatstr='%'+IDlen+'s%'+SeqIDlen+'s%'+Comp_IDlen+'s%'+Atom_IDlen\
                  +'s %'+Vallen+'.2f . 1\n'
    for linelist in data:
        ID=linelist[IDindex]
        SeqID=linelist[SeqIDindex]
        Comp_ID=linelist[Comp_IDindex]
        Atom_ID=linelist[Atom_IDindex]
        Val=float(linelist[Valindex])
        outfile.write(formatstr %(ID,SeqID,Comp_ID,Atom_ID,Val))
    infile.close()
    outfile.close()
    return outfilename


