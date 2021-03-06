import numpy as np


#Parametry modelu
popsize = 38386000
R_0 =  2.6     
sigma = 1/3.75 
gamma = 1/3.75 
omega1 = 1/8 
omega2 = 1/8 
epsilon1 = 0.05 
epsilon2 = 0.5 
epsilon3 = 0.5 
control  = 10 
kappa = 0.3 
seed = 30

#Początkowe warunki epidemii
inits = (popsize - 1,0,1,0,0,0,0,1)

#Raportowane dane
date =  np.array([4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
cum_cases = np.array([1,1,5,6,11,17,22,31,51,68,104,125,177,238,287,355,425,536,634,749,901,1051,1221])
cum_deaths = np.array([0,0,0,0,0,0,0,0,1,2,3,3,4,5,5,5,5,5,7,8,10,13,16])
inc_cases = np.array([1,0,4,1,5,6,5,9,20,17,36,21,52,61,49,68,70,111,98,115,152,150,170])
inc_deaths = np.array([0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,2,1,2,3,3])
no_tests = np.array([585,677,860,861,1157,1385,1652,2030,2238,2899,4425,5507,6699,7932,9556,11246,13119,15168,17678,20227,23025,26338,29700])


TEMP_R_0 = 2.6 
TEMP_TIME = None