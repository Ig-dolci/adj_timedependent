import numpy as np
import math as math
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=1000)


deltat = 1
timestep = 0.0001 
offset = 0
deltafd = 0.0001
tau = 4

#print(w)

arqp = '../base/DragLiftap.eny'
arqm = '../base/DragLiftam.eny'
    
    
Mp = np.loadtxt(arqp)

Mm = np.loadtxt(arqm)


grada = 0

init = 0
per = 240000

#print("FDm periodicity check", abs(Mm[init,1]-Mm[2000000,1]))
#print(Mm[init,1], Mm[2000000,1])

ft = 1;


for time in range(init, init+per):
    grada = grada + (Mp[time,1]*math.tanh(ft*(24-time*timestep*deltat))-Mm[time,1]*math.tanh(ft*(24-time*timestep*deltat))+Mp[time+1,1]*math.tanh(ft*(24-(time+1)*timestep*deltat))-Mm[time+1,1]*math.tanh(ft*(24-(time+1)*timestep*deltat)))/2*deltat*timestep/2/deltafd
    #grada = grada + (Mp[time,1]-Mm[time,1]+Mp[time+1,1]-Mm[time+1,1])/2*deltat*timestep/2/deltafd


  
    
#print('Newton-Cotes order 1 W')
print("grad FD A", grada)    

#grada = 0

#init = 2000000-40000
#per = 2000000-init

#for time in range(init, init+per,2):
#    grada = grada + (-Mm[time,1]-4*Mm[time+1,1]-Mm[time+2,1])/3*deltat*timestep/2/deltafd/(tau)
    
#init = 2000000-40000
#per = 2000000-init

#for time in range(init, init+per,2):
#    grada = grada + (Mp[time,1]+4*Mp[time+1,1]+Mp[time+2,1])/3*deltat*timestep/2/deltafd/(tau)
  
  
#print('Newton-Cotes order 2 W')
#print(grada)  





    




