import variables
import numpy as np
from scipy.integrate import odeint
import scipy

def model(x, t, seed, R_0, sigma, gamma, omega1, omega2, epsilon1, epsilon2, epsilon3, control, kappa):
    """Model SEIR
    S - podatni na zarażenie
    E - narażenia na kontakt
    I - zarażeni
    H - hospitalizowani
    V - poddani intensywnej terapii
    R - odporni po zetknięciu z chorobą
    D - zmarli
    C - całkowita liczba zarażonych
    """
    def beta(t,kappa,seed):
        """Parametr zmieniajcy się w momencie wprowadzenia rządowych regulacji"""
        if t < seed + variables.control:
            return R_0*gamma/variables.popsize
        else:
            return kappa*R_0*gamma/variables.popsize

    S,E,I,H,V,R,D,C = x 


    #---------- Implemetacja szumu --------------- 
    '''for time in variables.TEMP_TIME:
        if time == t:
            srand = np.random.normal(scale=I**1)
            I = I + srand
            S = S - srand
            C = C + srand'''

    dSdt = - beta(t,kappa,seed)*S*I         
    dEdt = beta(t,kappa,seed)*S*I - sigma*E
    dIdt = sigma*E - gamma*I 
    dHdt = epsilon1*gamma*I - omega1*H
    dVdt = epsilon2*omega1*H - omega2*V
    dRdt = (1 - epsilon1)*gamma*I + (1 - epsilon2)*omega1*H + (1 - epsilon3)*omega2*V
    dDdt = epsilon3*omega2*V
    dCdt = sigma*E
    return [dSdt ,dEdt ,dIdt ,dHdt ,dVdt ,dRdt ,dDdt ,dCdt]
#           0     1     2     3     4     5     6     7

def model_diff_solve(t, seed = variables.seed, kappa = variables.seed, R_0 = variables.R_0):
    """Funkcja rozwiązuje równanie różniczkowe dla modelu SEIR
    Argumenty:
    t - array
    tablica czasów, w których rozwiązuje się równania
    seed - float
    czas od pierwszego faktycznego zachorowania do do pierwszego odnotowanego przypadku
    kappa - float
    współczynik efektywności rządowych ograniczeń, znajduje się w przedziale od 0 do 1, gdzie 0 oznacza całkowite zatrzymanie rozwoju wirusa, a 1 brak wpływu ograniczeń
    """
    
    times = seed + t - t.min()
    times = np.append([0],times) # Roszerzamy tablicę czasów o 0, ponieważ dla tego czasu sa wartości variables.inits
    times = np.append(times , times.max()+1)
    variables.TEMP_TIME = times
    symulacja = odeint(model, variables.inits, times , args = (seed, R_0, variables.sigma, variables.gamma, variables.omega1, variables.omega2, variables.epsilon1, variables.epsilon2, variables.epsilon3, variables.control, kappa))
    return symulacja

def death_function(t,seed,kappa): 
    "Funkcja służąca do optymalizacji względem parametrów seed oraz kappa"
    times = seed + t - t.min()
    times = np.append([0],times) # Roszerzamy tablicę czasów o 0, ponieważ dla tego czasu sa wartości variables.inits
    times = np.append(times , times.max()+1)
    variables.TEMP_TIME = times
    symulacja = odeint(model, variables.inits, times , args = (seed, variables.TEMP_R_0, variables.sigma, variables.gamma, variables.omega1, variables.omega2, variables.epsilon1, variables.epsilon2, variables.epsilon3, variables.control, kappa))
    #---------------------------------------------------------------------------------------------------------------------------------
    #Część kodu odpowiedzialna, za wyłuskiwanie wartości iteresującej nas funkcji z tablicy rozwiązań symulacja, dla odpowednich chwil
    #---------------------------------------------------------------------------------------------------------------------------------
    D = []
    for i in symulacja[1:]:
        D.append(i[6])
    A = []
    for i in range(1,len(D)):
        A.append(D[i]-D[i-1])
    arrA = np.array(A)
    #---------------------------------------------------------------------------------------------------------------------------------
    return arrA


def data_fit(R_0 = 2.6):
    "Dopasowanie funkcji w zależności od parametru R_0"
    variables.TEMP_R_0 = R_0
    popt, pcov = scipy.optimize.curve_fit(death_function, variables.date, variables.inc_deaths, p0 = (variables.seed,variables.kappa), bounds = (0,[100,1])) # wartości parametrów oraz macierz kowariancji
    return popt, pcov