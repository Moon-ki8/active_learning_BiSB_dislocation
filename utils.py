import ase
import numpy as np
from ase.calculators.espresso import Espresso
from ase.io.espresso import write_espresso_in
import sys 
def data2xyz():
    after = ase.io.read('after.data',format='lammps-data',style='atomic')
    pos_after = after.positions
    natom_after = np.size(pos_after,0)
    species_after = after.get_chemical_symbols()
    species = [] 
    for sp in species_after:
        if sp == 'H':
            species.append('Bi')
        if sp == 'He':
            species.append('Sb')    
    after.set_chemical_symbols(species)
    fname = 'after.xyz'
    ase.io.write(fname, after, format ='extxyz')       

def eval_gap2xyz():
    quipoutfile = 'eval.xyz'
    # output .xyz file with gap error
    lineout = []
    with open('eval.out', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == 'A':
                line = line[3:]
                lineout.append(line)
        f.close()
    with open(quipoutfile, 'w') as f:
        f.writelines(lineout)
        f.close()
        
def find_max_gap_value():
    evalxyz = ase.io.read('EVAL/eval.xyz',format='extxyz')
    # Force error 
    f_var = abs(evalxyz.arrays['gap_variance_gradient'])
    f_std = np.sqrt(f_var)
    norm_std_ = np.linalg.norm(f_std,axis=1)
    #var = evalxyz.arrays['local_gap_variance']
    #value = np.max(var)
    pos = evalxyz.get_positions()
    # Find max value near dislocation core
    R = 100
    new_value,count,prev_v = [],0,0
    for curr_pos,value in zip(pos,norm_std_):
        dis = np.sqrt((curr_pos[0]-R)**2+(curr_pos[2]-R)**2)
        if dis < R-30 and value > prev_v:
            new_value.append(value)
            prev_v = value
            count += 1
            index = count
            pos_eval = curr_pos
    print(prev_v)
    return value 

def create_rep():
    file = 'relaxed.data'
    unit = ase.io.read(file,format='lammps-data',style='atomic')
    species_unit = unit.get_chemical_symbols()
    species = []
    for sp in species_unit:
        if sp == 'H':
            species.append('Bi')
        if sp == 'He':
            species.append('Sb')
    unit.set_chemical_symbols(species)
    # Import atom position and species 
    pos_unit = unit.positions
    natom_unit = np.size(pos_unit,0)
    species_unit = unit.get_chemical_symbols()
    # Compute axes 
    cell_unit = unit.cell
    x1_unit,x2_unit,x3_unit = cell_unit[0],cell_unit[1],cell_unit[2]
    # Import before
    file = 'before.data'
    bef = ase.io.read(file,format='lammps-data',style='atomic')
    species_bef = bef.get_chemical_symbols()
    pos_bef = bef.get_positions()
    cell_bef = bef.cell
    natom_bef = np.size(pos_bef,0)
    # Import cell ID 
    cellid = []
    with open('FBC_cellid.txt', 'r') as f:
        for line in f:
            temp = [int(item) for item in line.split()]
            cellid.append([int(temp[0]),int(temp[1]),int(temp[2])])
    # 1. Find highest var atom 
    file = 'eval.xyz'
    ori = ase.io.read(file,format='extxyz')
    atom_id = ori.arrays['id']
    grad = ori.arrays['local_gap_variance']
    #np.max(grad)
    #index_of_max = np.argmax(grad)
    pos = ori.get_positions()
    # Find max value near dislocation core
    R = 150 # !! THIS SHOULD BE CHANGED FOR SIZE 
    new_value,count,prev_v = [],0,0
    for curr_pos,value in zip(pos,grad):
        dis = np.sqrt((curr_pos[0]-R)**2+(curr_pos[2]-R)**2)
        if dis < R-8 and value > prev_v:
            new_value.append(value)
            prev_v = value
            count += 1
            index_of_max = count
    var_id = atom_id[index_of_max]
    curr_cell = cellid[var_id-1]
    #print(var_id)
    # In cell, find position  
    #     ------------
    #   / 1   /  2   /      z 
    #  / --- / ---  /      ^
    # /  4  /   3  /       | 
    # -------------        -> x
    pos_temp = pos_bef[var_id-1]- curr_cell[0]*x1_unit - curr_cell[2]*x3_unit
    
    zp = pos_temp[2]/x3_unit[2]
    pos_temp2  = pos_temp - zp*x3_unit
    yp = pos_temp2[1]/x2_unit[1]
    xp = pos_temp2[0]/x1_unit[0]
    if xp > 0.5 and zp > 0.5:
        flag = 2
    if xp < 0.5 and zp > 0.5:
        flag = 1
    if xp > 0.5 and zp < 0.5:
        flag = 3
    if xp < 0.5 and zp < 0.5:
        flag = 4
    
    if flag == 1:
        rep_cell_x,rep_cell_z = [-1,0],[0,1]
    if flag == 2:
        rep_cell_x,rep_cell_z = [0,1],[0,1]
    if flag == 3:
        rep_cell_x,rep_cell_z = [0,1],[-1,0]
    if flag == 4:
        rep_cell_x,rep_cell_z = [-1,0],[-1,0]
    # Find 2x2 atoms 
    cell_list = []
    for i in rep_cell_x:
        for j in rep_cell_z:
            cell_list.append([curr_cell[0]+i,0,curr_cell[2]+j])
    # Find atoms in cell_list
    local_id = []
    for i in cell_list:
        for count,j in enumerate(cellid):
            if i[0] == j[0] and i[2] == j[2]:
                #print(i)
                #print(j)
                local_id.append(count)
    # READ before and find
    file = 'after.data'
    curr = ase.io.read(file,format='lammps-data',style='atomic',sort_by_id=True)
    pos_curr = curr.get_positions()
    species_curr = curr.get_chemical_symbols()
    # Apply displacement and make periodic local image
    x1_trans = np.min([row[0] for row in cell_list])
    x2_trans = np.min([row[1] for row in cell_list])
    x3_trans = np.min([row[2] for row in cell_list])
    pos_rep,species_rep = [],[]
    #print(local_id)
    for i in local_id:
        dis0 = pos_curr[var_id-1] - pos_bef[var_id-1]
        dis = pos_curr[i] - pos_bef[i]
        distance = np.linalg.norm(pos_bef[i]-pos_bef[var_id-1])
        if distance < 4:
            position = pos_bef[i] + 1.0*dis - 1.0*dis0
        if distance < 10 and distance > 6:
            alpha = np.cos(np.pi/2/2*(distance-4))
            position = pos_bef[i] + 1.0*dis - 1.0*dis0
        else:
            position = pos_bef[i] + 1.0*dis - 1.0*dis0
        position += -x1_trans*x1_unit-x2_trans*x2_unit-x3_trans*x3_unit # Translation for cell
        pos_rep.append(position)
        if species_bef[i] == 'H':
            species_rep.append('Bi')
        if species_bef[i] == 'He':
            species_rep.append('Sb')
    cell_rep= cell_bef
    cell_rep[0],cell_rep[1],cell_rep[2] = x1_unit*2,x2_unit,x3_unit*2
    print(np.size(pos_rep))

    
    # Save as xyz file 
    rep = ase.io.read(fname,format='lammps-data',style='atomic')
    rep.set_chemical_symbols(species_rep)
    fname = 'rep.xyz'
    ase.io.write(fname, rep, format ='extxyz')

def DFT_rep():
    kpts = [4,2,2]
    ori= ase.io.read('rep.xyz',format='extxyz')
    pseudo_dict = {
       'Bi':'Bi.rel-pbe-dn-kjpaw_psl.1.0.0.UPF',
       'Sb':'Sb.rel-pbe-n-kjpaw_psl.1.0.0.UPF'
       }
    input_data = {
       'control':{
       'calculation':'scf',
       'pseudo_dir': '.',
       'tprnfor': True,
       'tstress': True,
       },
       'system':{
       'ecutwfc': 60,
       'ecutrho': 600,
       'degauss': 0.001,
       'occupations':'smearing',
       'smearing':'fd',
       'noncolin':True,
       'lspinorb':True,
       },
       'ELECTRONS':{
       'electron_maxstep':80,
       'conv_thr': 1e-6,
       'scf_must_converge':False,
       'mixing_beta':0.2
       }
    }
    calc = Espresso(pseudopotentials=pseudo_dict,input_data=input_data)
    ori.set_calculator(calc)
    # Write scf input file
    infile = open('scf.in', 'w')
    write_espresso_in(infile, ori, input_data=input_data, pseudopotentials=pseudo_dict,kpts=kpts)

def create_rep_new():
    file = 'unit_relaxed.data'
    unit = ase.io.read(file,format='lammps-data',style='atomic')
    species_unit = unit.get_chemical_symbols()
    species = []
    for sp in species_unit:
        if sp == 'H':
            species.append('Bi')
        if sp == 'He':
            species.append('Sb')
    unit.set_chemical_symbols(species)
    # Import atom position and species 
    pos_unit = unit.positions
    natom_unit = np.size(pos_unit,0)
    species_unit = unit.get_chemical_symbols()
    # Compute axes 
    cell_unit = unit.cell
    x1_unit,x2_unit,x3_unit = cell_unit[0],cell_unit[1],cell_unit[2]
    x1_norm,x2_norm,x3_norm = np.linalg.norm(x1_unit),np.linalg.norm(x2_unit),np.linalg.norm(x3_unit)

    # 1. Find highest var atom 
    file = 'eval.xyz'
    ori = ase.io.read(file,format='extxyz')
    atom_id = ori.arrays['id']
    #grad = ori.arrays['local_gap_variance']
    f_var = abs(ori.arrays['gap_variance_gradient'])
    f_std = np.sqrt(f_var)
    norm_std_ = np.linalg.norm(f_std,axis=1)
    pos = ori.get_positions()
    species = ori.get_chemical_symbols()
    # Find max value near dislocation core
    R = 100 # !! THIS SHOULD BE CHANGED FOR SIZE 
    new_value,prev_v = [],0
    for curr_pos,value,idd in zip(pos,norm_std_,atom_id):
        dis = np.sqrt((curr_pos[0]-R)**2+(curr_pos[2]-R)**2)
        if dis < R-30 and value > prev_v:
            new_value.append(value)
            prev_v = value
            index_of_max = idd
            
    index=int(np.where(atom_id == index_of_max)[0])
    #%% Create rep atom 
    rep_pos,rep_spe = [],[]
    #print(pos[index])
    for cpos,spe in zip(pos,species):
        rpos_ori = np.array([cpos[0]-pos[index][0],cpos[1],cpos[2]-pos[index][2]])
        #rpos_ori = cpos - pos[index] # Relative position
        n_rpos = np.linalg.norm(rpos_ori)
        if n_rpos == 0:
                rep_pos.append(rpos_ori)
                rep_spe.append(spe)
        if n_rpos > 0.0001 and n_rpos < 20 :
            zp = rpos_ori[2]/x3_unit[2]
            rpos  = rpos_ori - zp*x3_unit
            yp = rpos[1]/x2_unit[1]
            xp = rpos[0]/x1_unit[0]
            #print(xl,zl)
            if xp < 1 and xp > -1 and zp < 1 and zp > -1 :
                rep_pos.append(rpos_ori)
                rep_spe.append(spe)
    #print(rep_pos)
    rep_pos = rep_pos+x1_unit+x3_unit

    # Check overlap atom with images
    list_del = []
    for i in range(0,len(rep_pos)-1):
        for j in range(i+1,len(rep_pos)):
            cpos1,cpos2 = rep_pos[i],rep_pos[j]
            #print(cpos1-cpos2)
            for ii in [-1,0,1]:
                for jj in [-1,0,1]:
                    for kk in [-1,0,1]:
                        image_pos = cpos2+ii*2*x1_unit+jj*x2_unit+kk*2*x3_unit
                        if np.linalg.norm(cpos1-image_pos) < 1.5 and np.linalg.norm(cpos1-image_pos) > 0.001:
                            list_del.append([i])
                            #print(np.linalg.norm(cpos1-image_pos))
    # Your original list

    # Flatten the list and convert to a set to get unique elements
    unique_elements = set(elem for sublist in list_del for elem in sublist)

    # Convert the set back to a list
    list_del = list(unique_elements)
    if len(list_del)>0:
        rep_pos = np.delete(rep_pos, list_del, axis=0)
        rep_spe = np.delete(rep_spe, list_del, axis=0)

    rep = ase.Atoms(rep_spe,
                  positions=rep_pos,
                  cell=[2*x1_unit, x2_unit, 2*x3_unit],
                  pbc=[1, 1, 1])

    ase.io.write('rep.data',rep,format='lammps-data',atom_style='atomic')
    ase.io.write('rep.xyz',rep,format='extxyz')

#%% MAIN 
if len(sys.argv) ==  2:
        fun_name = sys.argv[1]
        
# Execute the function
#create_rep()
globals()[fun_name]()


