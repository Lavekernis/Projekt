import matplotlib.pyplot as plt
import model
import variables
import scipy
import numpy as np
from datetime import datetime, timedelta
from matplotlib import dates as mdates 

lista_R_0 = [2.0, 2.2, 2.4, 2.6]
kappa_list = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
global_view_kappa_list = [0,0.1,0.2,0.25,0.5,0.75,1]

def R_0_of_seed():
    fitted_seed = []
    error_seed = []
    for R_0 in lista_R_0:
       ppot ,pcov = model.data_fit(R_0 = R_0)
       fitted_seed.append(60 - ppot[0])
       error_seed.append(2*pcov[0,0]**0.5)

    plt.errorbar(fitted_seed ,lista_R_0 ,xerr= error_seed, fmt = 'rs')
    plt.xlabel("Data początku epidemii")
    plt.ylabel("Wartość $R_0$")

def R_0andKappa():
    for R_0 in lista_R_0:
        fig, ax = plt.subplots()
        ppot ,pcov = model.data_fit(R_0 = R_0)
        date_span = []
        date = datetime(2020,3,4)
        date_span.append(date)
        for i in range(len(variables.date)-1):
            date += timedelta(days = 1)
            date_span.append(date)
        
        for kappa in kappa_list:
            symulacja = model.model_diff_solve(variables.date, seed = ppot[0], kappa = kappa)
            I = []
            for i in symulacja[1:-1]:
                I.append(i[2])
            
            ax.plot_date(date_span,I, label = "Kappa " + str(kappa), fmt = '-')
            plt.xlabel("Czas")
            plt.ylabel("Liczba zarażonych")
            ax.legend()
        ax.plot_date(date_span,variables.cum_cases, fmt ='ok', label = "Raportowane przypadki")
        ax.legend()

def Global_View(kind,mouths):
    fig, ax = plt.subplots()
    ppot ,pcov = model.data_fit()
    time_span = np.linspace(int(ppot[0]),30*mouths,30*mouths)
    for kappa in global_view_kappa_list:
        symulacja = model.model_diff_solve(time_span,seed = ppot[0], kappa = kappa)
        if kind == 'H':
            a,label = 3, 'Hospitalizowani'
        elif kind == 'V':
            a,label = 4, 'Pacjenci na oddziale intensywnej terapii'
        elif  kind == 'D':
            a,label = 6, 'Zmarli'
        dane = []
        for i in symulacja[1:-1]:
            dane.append(i[a])
        date_span = []
        date = datetime(2020,3,4)
        date_span.append(date)
        for i in range(30*mouths-1):
            date += timedelta(days = 1)
            date_span.append(date)
        ax.plot_date(date_span, dane, label = "Kappa " + str(kappa), fmt = '-')
        plt.xlabel("Czas")
        plt.ylabel(label)
        ax.legend()

R_0_of_seed()
R_0andKappa()
#Global_View('D',10)
#Global_View('H',10)
#Global_View('V',10)
plt.show()