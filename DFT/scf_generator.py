import sys
import numpy as np
from ase.calculators.espresso import Espresso
from ase.io.espresso import write_espresso_in
import ase 
import ase.io
import random
from collections import Counter

# Specific seed number for reproducibility 
seed_number = 20
random.seed(seed_number)

def make_strained_cell_single(file,kpts):
    ori = ase.io.read(file,format='extxyz')
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
            'ecutwfc': 50,
            'ecutrho': 500,
            'degauss': 0.001,
            'occupations':'smearing',
            'smearing':'fd',
            'noncolin':True,
            'lspinorb':True,
            },
        'ELECTRONS':{
            'electron_maxstep':30,
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
    

kpts = [4,2,3]
file = 'rep.xyz'
make_strained_cell_single(file,kpts)

