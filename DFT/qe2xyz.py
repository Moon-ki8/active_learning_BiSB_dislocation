import ase 
import ase.io
temp_ATOM = ase.io.read('out.txt',index=0,format='espresso-out')
ase.io.write('rep_force.xyz',temp_ATOM,format='extxyz',append=False)

