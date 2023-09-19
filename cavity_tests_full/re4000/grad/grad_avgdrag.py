import numpy as np
import math as math
import subprocess
import sys
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=1000)


deltat = 1
timestep = 0.0001 
offset = 0
tau = 4

#print(w)

arqp = '../base/DragLift.eny'
    
Mp = np.loadtxt(arqp)


grada = 0

init = 0
per = 40000


print("Base flow periodicity check", abs(Mp[init,1]-Mp[init+per,1]))

print("J(t0) (term II)", Mp[0,1])

for time in range(init, init+per):
    grada = grada + (Mp[time,1]+Mp[time+1,1])/2*deltat*timestep
  
    
#print('Newton-Cotes order 1 W')
print("(term I)", grada)    

#grada = 0

#init = 0
#per = 40000

#for time in range(init, init+per,2):
#    grada = grada + (Mp[time,1]+4*Mp[time+1,1]+Mp[time+2,1])/3*deltat*timestep
  
  
#print('Newton-Cotes order 2 W')
#print("(term I)", grada)


grad_adj = subprocess.check_output(['tail', '-1', '../adjoint/adjgrad.fce']).decode(sys.stdout.encoding).strip().split()
print("Grad Adj A (term IV)", grad_adj[2])
print("Grad Adj B (term IV)", grad_adj[5])


    




