import os, shutil, random, string, subprocess

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

def seq2fasta(sequence,name='unknown',path='default',chain='A'):
    '''Generate the .fasta file needed
       sequence is a string of the residue sequence
       name is the name of the model
       path is where you want to save your .fasta file. The default is name.fasta
       you can specify the chain of the sequence. The default setting is chain A
    '''
    if path=='default':
        fasta_path=os.path.abspath(name+'.fasta')
    else:
        fasta_path=os.path.abspath(path)
    fastafile=open(fasta_path,'w')
    fastafile.write('>'+name+' | chain '+chain+'\n')
    for i in range(len(sequence)/60):
        fastafile.write(sequence[i*60:i*60+60]+'\n')
    if len(sequence)<60:
        i=-1
    fastafile.write(sequence[i*60+60:]+'\n')
    fastafile.close()
    return fasta_path

def make_fragments(sequence, name='unknown', fasta_path='default', chain='A', working_dir='default', cleanup=True, frag_3_path='default', frag_9_path='default', n_frags=200, n_candidates=1000, fragment_script='default', logger=False):
    ''' 
    Make fragments files for a given string sequence:
    The only necessary parameter is sequence but passing in the name is also encouraged as you don't want
    a structure to be called 'unknown'
    Before using this utility. Make sure you have everything set correctly in the make_fragment script in
    rosetta_tools.
    
    sequence: the squence of the protein
    name: the name of the protein. 
    fasta_path: the path where you want to save your fasta file. By default, the fasta file will save to
                your current directory and will have the name of the protein with the extension .fasta.
                The path is a path to the file not to the folder containing the file.
    chain: the chain of the sequence. Default is A
    working_dir: The working dir where temporary folders will be creatd. Default is './working'. If this
                 folder does not exist, it will be created.
    cleanup: whether to delete the temporary folders or not. The default is True.
    frag_3_path: The file path to the file where you want to save your 3 piece fragments. By default, it
                 will be saved to current folder with the name of your protein and an extension of .3mers
                 Usually, you will not worry too much about the path as it is passed out by function return.
    frag_9_path: The file path to the file where you want to save your 9 piece fragments. Having the same
                 naming rule with frag_3_path.
    n_frags: number of framgents to be generated per residue. Default setting is 200
    n_candidates: number of fragment candidates. Default setting is 1000
    fragment_script: the location of the script. By default it is $COMPASSInstall/rosetta3.4/rosetta_tools/
                     make_fragments.pl

    Function returns:
    (fasta_file, frag_3_file, frag_9_file)
    fasta_file: a string path to the fasta file
    frag_3_file: a string path to the 3 piece fragments file
    frag_9_file: a string path to the 9 piece fragments file    
    ''' 
    #mark the initial folder. The program will change dir to the working dir but will cd back at the end of the function
    initial_dir=os.getcwd()
    
    fasta_file=seq2fasta(sequence, name=name, path=fasta_path, chain=chain)
    if frag_3_path=='default':
        frag_3_file=os.path.abspath('%s.3mers' %name)
    else:
        frag_3_file=os.path.abspath(frag_3_path)

    if frag_9_path=='default':
        frag_9_file=os.path.abspath('%s.9mers' %name)
    else:
        frag_9_file=os.path.abspath(frag_9_path)
    
    if fragment_script=='default':
        script_path=os.environ['COMPASSInstall']+'/rosetta3.4/rosetta_tools/fragment_tools/make_fragments.pl'
    else:
        script_path=os.path.abspath(fragment_script)
         
    #check/create working_dir and save abs path to working_path
    if working_dir=='default':
        working_path='./working'
    else:
        working_path=working_dir
    if not os.path.exists(working_path):
        os.mkdir(working_path)
    working_path=os.path.abspath(working_path)+'/'
    
    #generate the temporary folder name
    while True:
        temp_ID=''.join(random.choice(string.ascii_uppercase) for x in range(4))+chain
        temp_dir=working_path+temp_ID
        if not os.path.exists(temp_dir):
            break
    os.mkdir(temp_dir)
    os.chdir(temp_dir)
    argument=[script_path, '-id', temp_ID, '-n_frags', str(n_frags), '-n_candidates',str(n_candidates),fasta_file]
    
    if logger:
        frag_process = subprocess.Popen(args=argument, stderr = subprocess.STDOUT, stdout = subprocess.PIPE)
    else:
        frag_process = subprocess.Popen(args=argument)

    print 'Process %d command: ' % frag_process.pid,  ' '.join(argument)

    
    frag_process.wait()
    if logger:
        print "####### output for process %d #######" % frag_process.pid
        print frag_process.stdout.read()
    
    temp_3_frag=os.path.abspath('%s.%d.3mers' %(temp_ID,n_frags))
    temp_9_frag=os.path.abspath('%s.%d.9mers' %(temp_ID,n_frags))

    os.rename(temp_3_frag, frag_3_file)
    os.rename(temp_9_frag, frag_9_file)
    
    os.chdir(initial_dir)

    if cleanup:
        shutil.rmtree(temp_dir)
    return(fasta_file,frag_3_file,frag_9_file) 

def generate_models(name='default', fasta_path='default', frag_3_path='default', frag_9_path='default', model_dir='default', working_dir='default', num_model=100, num_process=8, pdb_name_list=[], rosetta_abinitio='default', rosetta_DB='default', logger = False):
    '''
    Generate Structure from given fasta, frag_3 and frag_9 files. 
    You can run the program by specify the name and the program will search for the needed
    files with the name in the current dir. Or you can specify frag_3_path and frag_9_path and 
    the program will load those two files. When both specified, function will give preference
    to what frag_3_path and frag_9_path

    name: name of the protein
    fasta_path: file path to the fasta file
    frag_3_path: file path to the 3 pieces fragments file
    frag_9_path: file path to the 9 pieces fragments file
    model_dir: folder path to the dir where you want to save your models. Default is ./models
    working_dir:folder  path to the working dir where the function will create temporary files. Default is ./working
    num_model: num of models you want to generate. Default is 100
    num_process: num of process to be created. Default is 8
    pdb_name_list: specifing the name of the models. If you specify the list, it must have the length of
                   num_model. If not specified the models will have the name rosetta_model_00000001.pdb
                   and such.
    rosetta_abinitio: the abinitio program you will run. The default is
                      $COMPASSInstall/rosetta3.4/rosetta_source/bin/AbinitioRelax.macosgccrelease
    rosetta_DB: the rosetta Database. THe default is
                $COMPASSInstall/rosetta3.4/rosetta_database/

    Function return:
    A list of paths to the generated pdb files.
    '''
    if name!='default':
        fasta_file=os.path.abspath(name+'.fasta')
        frag_3_file=os.path.abspath(name+'.3mers')
        frag_9_file=os.path.abspath(name+'.9mers')
    if fasta_path!='default':
        fasta_file=os.path.abspath(fasta_path)
    if frag_3_path!='default':
        frag_3_file=os.path.abspath(frag_3_path)
    if frag_9_path!='default':
        frag_9_file=os.path.abspath(frag_9_path)
    
    if not os.path.exists(fasta_file):
        raise(Exception('Fasta file not found'))
    if not os.path.exists(frag_3_file):
        raise(Exception('3 piece fragments file not found'))
    if not os.path.exists(frag_9_file):
        raise(Exception('9 piece fragments file not found'))

    if model_dir=='default':
        model_path=os.path.abspath('./models')
    else:
        model_path=os.path.abspath(model_dir)
    if not os.path.exists(model_path):
        os.mkdir(model_path)

    if working_dir=='default':
        working_path=os.path.abspath('./working')
    else:
        working_path=os.path.abspath(working_dir)
    if not os.path.exists(working_path):
        os.mkdir(working_path)

    if rosetta_abinitio=='default':
        abinitio_path=os.environ['COMPASSInstall']+'/rosetta3.4/rosetta_source/bin/AbinitioRelax.macosgccrelease'
    else:
        abinitio_path=os.path.abspath(rosetta_abinitio)
    if not os.path.exists(abinitio_path):
        raise(Exception('Abinitio Relax program not found!'))

    if rosetta_DB=='default':
        database_dir=os.environ['COMPASSInstall']+'/rosetta3.4/rosetta_database/'
    else:
        database_dir=os.path.abspath(rosetta_DB)
    if not os.path.exists(database_dir):
        raise(Exception('Rosetta Database not found!'))
    if pdb_name_list==[]:
        model_file_list=[model_path+'/rosetta_model_%s.pdb' %str(x+1).zfill(8) for x in range(num_model)]
    else:
        if len(pdb_name_list)!=num_model:
            raise(Exception('Specified pdb name list length does not math number of models'))
        model_file_list=[model_path+'/'+pdbname for pdbname in pdb_name_list]
    
    #divide jobs
    #processlist is a list of Job objects
    initial_path=os.getcwd()
    processlist=[]
    for i in range(num_process):
        processlist.append(Job())
    counter=0
    for model in model_file_list:
        if counter==num_process:
            counter=0
        processlist[counter].models.append(model)
        counter+=1

    #in case you have less than num_process of models 
    for i in range(len(processlist)):
        if processlist[i].models==[]:
            processlist=processlist[:i]
            break
    #creating temp folders for the process
    for i in range(len(processlist)):
        while True:
            temp_folder='model_temp_'+''.join(random.choice(string.ascii_uppercase) for x in range(6))
            temp_folder=working_path+'/'+temp_folder
            if not os.path.exists(temp_folder):
                break
        os.mkdir(temp_folder)
        processlist[i].folder=(temp_folder)
    for i in range(len(processlist)):
        os.chdir(processlist[i].folder)
        argument=[abinitio_path,'-in::file::fasta',fasta_file,'-in:file:frag3',frag_3_file,'-in:file:frag9',frag_9_file,\
                  '-database',database_dir,'-nstruct',str(len(processlist[i].models)),'-out:pdb', '-out:path',processlist[i].folder,\
                  '-constant_seed','-jran',str(random.randint(1000000,9999999))]


        if logger:
            singleprocess = subprocess.Popen(args=argument, stderr = subprocess.STDOUT, stdout = subprocess.PIPE)
        else:
            singleprocess = subprocess.Popen(args=argument)
        processlist[i].popen=singleprocess
    os.chdir(initial_path)
    for i in range(len(processlist)):
        processlist[i].popen.wait()
    for i in range(len(processlist)):
        if logger:
            print "####### output for process %d #######" % processlist[i].pid
            print processlist[i].stdout.read()
        for j in range(len(processlist[i].models)):
            temp_pdb=processlist[i].folder+'/S_%s.pdb'%str(j+1).zfill(8)
            os.rename(temp_pdb,processlist[i].models[j])
        shutil.rmtree(processlist[i].folder)
    return model_file_list

def runrosetta(sequence,name='unknown',fasta_path='default',chain='A',working_dir='default',model_dir='default',cleanup=True, frag_3_path='default', frag_9_path='default',n_frags=200, n_candidates=1000, num_model=128, num_process=8, pdb_name_list=[], rosetta_abinitio='default', rosetta_DB='default', fragment_script='default'):
    '''
    This function calls make_fragment and generate_models to create all the models you need.
    You must provide the sequence and providing the name is recommended. Providing the file saving location
    is option. By default generate 128 models with 8 processes.

    sequence: the squence of the protein
    name: the name of the protein. 
    fasta_path: the path where you want to save your fasta file. By default, the fasta file will save to
                your current directory and will have the name of the protein with the extension .fasta.
                The path is a path to the file not to the folder containing the file.
    chain: the chain of the sequence. Default is A
    working_dir: The working dir where temporary folders will be creatd. Default is './working'. If this
                 folder does not exist, it will be created.
    cleanup: whether to delete the temporary folders or not. The default is True.
    frag_3_path: The file path to the file where you want to save your 3 piece fragments. By default, it
                 will be saved to current folder with the name of your protein and an extension of .3mers
                 Usually, you will not worry too much about the path as it is passed out by function return.
    frag_9_path: The file path to the file where you want to save your 9 piece fragments. Having the same
                 naming rule with frag_3_path.
    n_frags: number of framgents to be generated per residue. Default setting is 200
    n_candidates: number of fragment candidates. Default setting is 1000
    fragment_script: the location of the script. By default it is $COMPASSInstall/rosetta3.4/rosetta_tools/
                     make_fragments.pl
    model_dir: folder path to the dir where you want to save your models. Default is ./models
    num_model: num of models you want to generate. Default is 100
    num_process: num of process to be created. Default is 8
    pdb_name_list: specifing the name of the models. If you specify the list, it must have the length of
                   num_model. If not specified the models will have the name rosetta_model_00000001.pdb
                   and such.
    rosetta_abinitio: the abinitio program you will run. The default is
                      $COMPASSInstall/rosetta3.4/rosetta_source/bin/AbinitioRelax.macosgccrelease
    rosetta_DB: the rosetta Database. THe default is
                $COMPASSInstall/rosetta3.4/rosetta_database/


    Function returns:
    (fasta_file, frag_3_file, frag_9_file, model_paths)
    fasta_file: a string path to the fasta file
    frag_3_file: a string path to the 3 piece fragments file
    frag_9_file: a string path to the 9 piece fragments file    
    model_paths: a list of paths to the models
    '''
    (fasta_file, frag3_file,frag9_file) = make_fragment(sequence=sequence, name=name, fasta_path=fasta_path,\
                                                       chain=chain, working_dir=working_dir,cleanup=cleanup,\
                                                       frag_3_path=frag_3_path, frag_9_path=frag_9_path,\
                                                       n_frags=n_frags, n_candidates=n_candidates, \
                                                       fragment_script=fragment_script)
    model_list=generate_models(fasta_path=fasta_file, frag_3_path=frag3_file, frag_9_path=frag9_file, model_dir=model_dir,\
                               working_dir=working_dir, num_model=num_model, num_process=num_process, pdb_name_list=pdb_name_list,\
                               rosetta_abinitio=rosetta_abinitio, rosetta_DB=rosetta_DB)
    return (fasta_file, frag3_file, frag9_file, model_list)

def relax1model(model,output='default',relaxprogram='default',rosettaDB='default',WAIT=True,fast=False, logger = None ):
    ''' 
    relax one pdb model
    
    parameters:
    model : the path to one model you want to relax.
    output : the path to the output file. By default, this will be adding '_relaxed' to the end of the pdb file name.
            This option is only useful when WAIT is True
    relaxprogram : the path to the relax program. The default is
                   $COMPASSInstall/rosetta3.4/rosetta_source/bin/relax.macosgccrelease
    rosettadb : the path to the Rosetta database. The default is
                $COMPASSInstall/rosetta3.4/rosetta_database
    WAIT: a bool sepcifing if you want the function to wait the relaxiation process to finish. The default is True.
          In the wait mode, the output file is also renamed correctly to your intended output. Otherwise, it will be
          renamed to the output_path returned by the function.
    fast: a bool indicating to do a fast or a thorough relax. The default if false.
    Returns:
    (output_path, POPEN)
    output_path is a string. It stores the path to the output file
    POPEN is the subprocess.Popen object used to control the process
    '''
    
    model_path=os.path.abspath(model)

    #setting rosetta relax and rosetta database
    if relaxprogram=='default':
        program_path=os.environ['COMPASSInstall']+'/rosetta3.4/rosetta_source/bin/relax.macosgccrelease'
    else:
        program_path=os.path.abspath(relaxprogram)
    
    if rosettaDB=='default':
        database_path=os.environ['COMPASSInstall']+'/rosetta3.4/rosetta_database'
    else:
        database_path=os.path.abspath(rosettaDB)
    
    if fast:
        fastflag='-relax:fast'
    else:
        fastflag='-relax:thorough'

    #creating subprocess
    argument=[program_path,'-database',database_path,'-in:file:s',model_path,'--nstruct','1',fastflag]
    if logger:
        POPEN = subprocess.Popen(args=argument, stderr = subprocess.STDOUT, stdout = subprocess.PIPE)
    else:
        POPEN = subprocess.Popen(args=argument)

    
    if output=='default':
        if model_path.endswith('.pdb'):
            output_path=model_path[:-4]+'_relaxed.pdb'
        else:
            output_path=model_path+'_relaxed.pdb'
    else:
        output_path=os.path.abspath(output)
   
    if WAIT:
        POPEN.wait()
        if logger:
            print "####### output for process %d #######" % POPEN.pid
            print POPEN.stdout.read()
        os.rename(os.path.basename(model_path)[:-4]+'_0001.pdb',output_path)
    else:
        output_path=os.path.abspath(os.path.basename(model_path)[:-4]+'_0001.pdb')

    return (output_path,POPEN)

def relaxlist(model_list,out_list=[],relaxprogram='default',rosettaDB='default',fast=False,CPUnum=8, logger = False):
    '''
    This function will relax all the pdb files in the model_list
    model_list : a list of strings of paths to the input pdb files.
    out_list : optional. When specified, it is the list of stings of paths to the output files.
    If you want to sepcify output_list, make sure that it has the same length with model_list
    relaxprogram : the path to the relax program. The default is
                   $COMPASSInstall/rosetta3.4/rosetta_source/bin/relax.macosgccrelease
    rosettadb : the path to the Rosetta database. The default is
                $COMPASSInstall/rosetta3.4/rosetta_database
    fast: a bool indicating to do a fast or a thorough relax. The default if false.

    Function return:
    The function returns a list of paths to the output file
    '''
    output_list=[]
    if out_list==[]:
        for path in model_list:
            if path.endswith('.pdb'):
                output_list.append(os.path.abspath(path[:-4]+'_relaxed.pdb'))
            else:
                output_list.append(os.path.abspath(path+'_relaxed.pdb'))
    else:
        if not len(model_list)==len(output_list):
            print("ERROR: Length of output_list must match length of model_list")
        output_list=[os.path.abspath(path) for path in out_list]
    temp_path_list=[]
    POPEN_list=[]
    for i in range(len(model_list)):
        (temp_path,POPEN)=relax1model(model_list[i],relaxprogram=relaxprogram,rosettaDB=rosettaDB,fast=fast,WAIT=False,logger=logger)
        temp_path_list.append(temp_path)
        POPEN_list.append(POPEN)
        print "Popen object created with PID:"
        print POPEN.pid
        if ((i+1)%CPUnum==0) or (i==len(model_list)-1):
            for j in range(max(0,i-CPUnum+1),i+1):
                POPEN_list[j].wait()
            if logger:
                print "####### output for process %d #######" % POPEN_list[j].pid
                print POPEN_list[j].stdout.read()

    for i in range(len(model_list)):
        os.rename(temp_path_list[i],output_list[i])

    return output_list
