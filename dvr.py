import numpy as np
import subprocess
from numpy import pi
from numpy import *
import matplotlib.pyplot as plt

dist_range = np.arange(0.6572,1.7572,0.04)

h_bar = 1 # A.U 

delta_x =  0.07558904535685 # A.U

cf = (1.66054e-27/9.10938356e-31)
#m = 1.007825*cf # u Zeroth order(mass of H) 

m = 0.9514424874*cf #1st order (reduced mass of OH and H) 
 
const = ((h_bar)**2)/((2*m)*(delta_x)**2)

penergy = []
energy_matrix = []

mu_x = []
mu_y = []
mu_z = []
mu = []


for i,dist in enumerate(dist_range):

      filename = "new" + str(i) + ".log"     
      s = "grep 'SCF Done:' %s "%(filename)
      out = subprocess.check_output(s, shell=True)
      line = out.split()
      pot_energy = float(line[4])
      penergy.append(pot_energy)
      column =[]
      for j,dist in enumerate(dist_range):
         if i==j:
           kin_energy = const*((-1)**(i-j))*((pi**2)/3)
           E = pot_energy + kin_energy
           column.append(E)
         else:
              kin_energy = const*((-1)**(i-j))*(2/((i-j)**2))
              column.append(kin_energy)
      energy_matrix.append(column)
      
      t = "grep 'Dipole moment (field-independent basis, Debye):' %s -A 1"%(filename)
      out1 = subprocess.check_output(t, shell=True)
      line1 = out1.split()
      u_x = float(line1[6])*0.393430307
      u_y = float(line1[8])*0.393430307
      u_z = float(line1[10])*0.393430307
      u = float(line1[12])*0.393430307
      
      mu_x.append(u_x)
      mu_y.append(u_y)
      mu_z.append(u_z)
      mu.append(u)



energy_ev, energy_es = np.linalg.eigh(energy_matrix)

print ("Eigenenergies for 20  molecules in the cluster ",energy_ev)
print (mu_x,mu_y,mu_z,mu) 
print ("states" , energy_es)

freq = energy_ev[1]-energy_ev[0] #A.U.
freq_cm = freq*219474.6305 #cm-1


print ("Transition frequency",freq_cm)
one = energy_es[:,1]
zero = energy_es[:,0]

print ("Eigenstate(0) = ", zero)
print ("Eigenstate(1) = " , one)

index = np.argmin(penergy)
print (index)

req = dist_range[index]*1.8897259886
print(req)

xmatrix = np.zeros((28,28))

for k, dist in enumerate(dist_range):

#    c = []
    for n, d in enumerate(dist_range):
       if k==n:
         value = req - d
         xmatrix[k,n]= value*1.8897259886
       else:
           continue
#       c.append(value)
#    xmatrix.append(c)

print ('xmatrix = ',xmatrix) 


x10 = np.linalg.multi_dot([zero,xmatrix,one])
print ('x10 = ',x10) 

#du_x = (mu_x[index-2] - 8.0*mu_x[index-1] + 8.0*mu_x[index+1] - mu_x[index+2])/ (12.0*delta_x)
#du_y = (mu_y[index-2] - 8.0*mu_y[index-1] + 8.0*mu_y[index+1] - mu_y[index+2])/ (12.0*delta_x)
#du_z = (mu_z[index-2] - 8.0*mu_z[index-1] + 8.0*mu_z[index+1] - mu_z[index+2])/ (12.0*delta_x)



du_x = (mu_x[index+1] - mu_x[index-1])/(2.0*delta_x)
du_y = (mu_y[index+1] - mu_y[index-1])/(2.0*delta_x)
du_z = (mu_z[index+1] - mu_z[index-1])/(2.0*delta_x) 
du = np.sqrt( du_x**2 + du_y**2 + du_z**2)    

print ("du =  ", du)

td = du*x10
print ("Transition dipole = ", td)


