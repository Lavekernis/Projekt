

#parms ={'seed':, 'R_0':, 'sigma':, 'gamma':, 'omega1':, 'omega2':, 'epsilon1':, 'epsilon2':, 'epsilon3':, 'control':, 'kappa':}

#fixed = {'seed':  np.log(30), #seed temp
         # 'R_0':  2.2,           #'R_0'
        #   'sigma': 1/3.75, #'sigma'
       #     'gamma': 1/3.75, #'gamma'
      #      'omega1': 1/8, #'omega1'
     #       'omega2': 1/8, #'omega2'
    #       'epsilon1': 0.05, # 'epsilon1'
   #        'epsilon2': 0.5, # 'epsilon2'
  #         'epsilon3': 0.5, # 'epsilon3'
 #          'control': 10, #  'control'
#'kappa': scipy.special.logit(0.1)} # kappa temp

#    0     1   2       3       4       5       6           7           8       9       10
#seed, R_0, sigma, gamma, omega1, omega2, epsilon1, epsilon2, epsilon3, control, kappa

#def trans(pars):
 #   pars[0] = math.exp(pars[0])
  #  pars[10] = scipy.special.expit(pars[10])
   # return tuple(pars)

#def nll(seed, R_0, sigma, gamma, omega1, omega2, epsilon1, epsilon2, epsilon3, control, kappa):
 #   param = [seed, R_0, sigma, gamma, omega1, omega2, epsilon1, epsilon2, epsilon3, control, kappa]
  #  param = trans(param)
   # times = seed + date - date.min()
    #times = np.append([0],times)
    #times = np.append(times, times.max()+1)
    #symulacja = odeint(model, inits, [1,2], args = param)
    #symulacja = symulacja[1:]
    #return symulacja



import matplotlib.pyplot as plt
import math
import numpy as np
import scipy
from scipy.integrate import odeint

popsize = 38386000

R_0 =  2.6      #'R_0'
sigma = 1/3.75 #'sigma'
gamma = 1/3.75 #'gamma'
omega1 = 1/8 #'omega1'
omega2 = 1/8 #'omega2'
epsilon1 = 0.05 # 'epsilon1'
epsilon2 = 0.5 # 'epsilon2'
epsilon3 = 0.5 # 'epsilon3'
control  = 10 #  'control'

inits = (popsize - 1,0,1,0,0,0,0,1)

date =  np.array([4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
cum_cases = np.array([1,1,5,6,11,17,22,31,51,68,104,125,177,238,287,355,425,536,634,749,901,1051,1221])
cum_deaths = np.array([0,0,0,0,0,0,0,0,1,2,3,3,4,5,5,5,5,5,7,8,10,13,16])
inc_cases = np.array([1,0,4,1,5,6,5,9,20,17,36,21,52,61,49,68,70,111,98,115,152,150,170])
inc_deaths = np.array([0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,2,1,2,3,3])
no_tests = np.array([585,677,860,861,1157,1385,1652,2030,2238,2899,4425,5507,6699,7932,9556,11246,13119,15168,17678,20227,23025,26338,29700])







def model(x, t, seed, R_0, sigma, gamma, omega1, omega2, epsilon1, epsilon2, epsilon3, control, kappa):
    
    def beta(t,kappa,seed):
        if t < seed + control:
            return R_0*gamma/popsize
        else:
            return kappa*R_0*gamma/popsize

    S,E,I,H,V,R,D,C = x
    dS = - beta(t,kappa,seed)*S*I
    dE = beta(t,kappa,seed)*S*I - sigma*E
    dI = sigma*E - gamma*I 
    dH = epsilon1*gamma*I - omega1*H
    dV = epsilon2*omega1*H - omega2*V
    dR = (1 - epsilon1)*gamma*I + (1 - epsilon2)*omega1*H + (1 - epsilon3)*omega2*V
    dD = epsilon3*omega2*V
    dC = sigma*E
    return [dS,dE,dI,dH,dV,dR,dD,dC]
#           0  1  2  3   4  5  6  7



def death_function(t,seed,kappa): 
    times = seed + t - t.min()
    times = np.append([0],times)
    times = np.append(times , times.max()+1)
    symulacja = odeint(model, inits, times , args = (seed, R_0, sigma, gamma, omega1, omega2, epsilon1, epsilon2, epsilon3, control, kappa))
    D = []
    for i in symulacja[1:]:
        D.append(i[6])
    A = []
    for i in range(1,len(D)):
        A.append(D[i]-D[i-1])
    arrA = np.array(A)
    return arrA




#np.log(30) #seed temp
kappa = 0.3 #scipy.special.logit(0.1) # kappa temp
seed = 30
popt, pcov = scipy.optimize.curve_fit(death_function, date, inc_deaths, p0 = (seed,kappa), bounds = (0,[100,1])) 
print(popt)
print(pcov)
plt.show()