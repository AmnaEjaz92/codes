import numpy as np
import mdtraj as md
import os
import subprocess

def minImage(d, ucl):
   return np.abs(d - ucl*round(d/ucl))

def minImage_(d, ucl):
   return (d - ucl*round(d/ucl))

def norm(v):
   return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

# run 'vpkg_require python/20200808:spec' to load mdtraj

gromacs_file = '/work/akanane/users/Amna/acn_water/trajectory/prod.gro'
xtc_file = '/work/akanane/users/Amna/acn_water/trajectory/prod_nopbc.xtc'

# Loading entire trajectory
traj = md.load(xtc_file, top=gromacs_file)

xyz = 10.0*traj.xyz
ucl = 10.0*traj.unitcell_lengths
nfr = traj.n_frames
natoms = traj.n_atoms
topology = traj.topology

print ("unit cell lengths",ucl[0,:])
print ("The center of the box is at",ucl[0,:]*(0.5))

# ACN molecule coordinates and Center of mass of CN

c_xyz = xyz[0,4,:]
n_xyz = xyz[0,5,:]
mN = 14.0067
mC = 12.0107

com_CN = np.zeros((3))
com_CN[0] = ((mN*n_xyz[0]+mC*c_xyz[0])/(mN+mC))
com_CN[1] = ((mN*n_xyz[1]+mC*c_xyz[1])/(mN+mC))
com_CN[2] = ((mN*n_xyz[2]+mC*c_xyz[2])/(mN+mC))


print ("Coordinates of Carbon", c_xyz)
print ("Coordinates of Nitrogen", n_xyz)
print ("Coordinates of COM", com_CN)

#Unit vector along the stretch
v = n_xyz - c_xyz
d = norm(v)
unit_vector = v/d

#Grid
grid  = np.arange(0.5572,1.8572,0.02)
l = len(grid)
print ("Length of Grid", l)

# Water molecules
water_molecules = [6,2005,4]

#Sorting data
data = []
for molecule in water_molecules:

      o_xyz_loop = xyz[0,molecule,:]

      vv = o_xyz_loop - com_CN

      r_x = minImage_(vv[0], ucl[0,0])
      r_y = minImage_(vv[1], ucl[0,1])
      r_z = minImage_(vv[2], ucl[0,2])

      cd = np.sqrt(r_x**2 + r_y**2 + r_z**2)

      loop_data = []
      loop_data.append(water_molecule)
      loop_data.append(cd)

      data.append(loop_data)   

data_sorted = sorted(data, key=lambda x: x[1])

#Choosing cluster around Nitrogen

closest = []
for i in range(8):
   name1 = data_sorted[i]
   mol_index = name1[0]
   distance = name1[1]
   closest.append(mol_index)

#Point charges

point_charges = []
for k in range(8,68):
   name2 = data_sorted[k]
   mol_indx = name2[0]
   distance2 = name2[1]
   point_charges.append(mol_indx)

#Gaussian Input Files
   
for grid_id, dist in enumerate(grid):
      filename = str(grid_id) + ".com" 
      n_new_xyz = c_xyz + dist*unit_vector      
      f = open(filename, 'w')
      f.write('# B3LYP/6-311++G(d,p) NoSymm Int=Ultrafine Charge SCF=tight Test \n')
      f.write('\n')
      f.write('comment line\n')
      f.write('\n')
      f.write('0 1\n')
      f.write('%s  %7.5f  %7.5f   %7.5f \n'%("C",xyz[0,0,0],xyz[0,0,1],xyz[0,0,2]))
      f.write('%s  %7.5f  %7.5f   %7.5f \n'%("H",xyz[0,1,0],xyz[0,1,1],xyz[0,1,2]))
      f.write('%s  %7.5f  %7.5f   %7.5f \n'%("H",xyz[0,2,0],xyz[0,2,1],xyz[0,2,2]))
      f.write('%s  %7.5f  %7.5f   %7.5f \n'%("H",xyz[0,3,0],xyz[0,3,1],xyz[0,3,2]))
      f.write('%s  %7.5f  %7.5f   %7.5f \n'%("C",c_xyz[0],c_xyz[1],c_xyz[2]))
      f.write('%s  %7.5f  %7.5f   %7.5f \n'%("N",n_new_xyz[0],n_new_xyz[1],n_new_xyz[2]))
      for j in closest:
         e = j+4
         xyz_20 = xyz[0,j:e,:]
         o_20 = xyz_20[0,:]
         h1_20 = xyz_20[1,:]
         h2_20 = xyz_20[2,:]
   # correct coordinates here
         voo = np.zeros((3))     
         voo[:] = com_CN[:] - o_20[:]
         sv = np.zeros((3))     
         sv[0] = (voo[0] - minImage_(voo[0], ucl[0,0]))
         sv[1] = (voo[1] - minImage_(voo[1], ucl[0,1]))
         sv[2] = (voo[2] - minImage_(voo[2], ucl[0,2]))
 
         f.write('%s  %7.5f  %7.5f   %7.5f \n'%("O",o_20[0]+sv[0],  o_20[1]+sv[1],  o_20[2]+sv[2]))
         f.write('%s  %7.5f  %7.5f   %7.5f \n'%("H",h1_20[0]+sv[0], h1_20[1]+sv[1], h1_20[2]+sv[2]))
         f.write('%s  %7.5f  %7.5f   %7.5f \n'%("H",h2_20[0]+sv[0], h2_20[1]+sv[1], h2_20[2]+sv[2]))

      f.write('\n')
      for m in point_charges:   
         ec = m + 4
         xyz_c = xyz[0,m:ec,:]
         h1c = xyz_c[1,:]
         h2c = xyz_c[2,:]
         mc = xyz_c[3,:]
         oc = xyz_c[0,:]
         vooc = np.zeros((3))     
         vooc[:] = oc[:] - com_CN[:]

         svc = np.zeros((3))     
 
         svc[0] = (vooc[0] - minImage_(vooc[0], ucl[0,0]))
         svc[1] = (vooc[1] - minImage_(vooc[1], ucl[0,1]))
         svc[2] = (vooc[2] - minImage_(vooc[2], ucl[0,2]))
 

         f.write('%7.5f  %7.5f  %7.5f  %3.2f \n'%(h1c[0]+svc[0],h1c[1]+svc[1],h1c[2]+svc[2],0.52))
         f.write('%7.5f  %7.5f  %7.5f  %3.2f \n'%(h2c[0]+svc[0],h2c[1]+svc[1],h2c[2]+svc[2],0.52))
         f.write('%7.5f  %7.5f  %7.5f  %3.2f \n'%(mc[0]+svc[0],mc[1]+svc[1],mc[2]+svc[2],-1.04))
      f.write('\n')
      f.close() 
   


