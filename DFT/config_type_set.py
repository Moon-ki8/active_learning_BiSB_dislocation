import ase.io
temp=ase.io.read('rep_force.xyz',format='extxyz')
temp.info['config_type']='rep'
ase.io.write('rep_force.xyz',temp,format='extxyz')
