import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import log as ln
from numpy import exp
import os
import math
from datetime import datetime  
import time

hi_values = pd.read_csv('./hi.csv', sep = ';').to_numpy() #importa os valores de hi do arquivo csv e transforma em uma matriz 
bonding_alpha = pd.read_csv('./bonding_alpha.csv', sep = ';').to_numpy() #importa os valores da energia de ligação da estrutura alpha
bonding_delta = pd.read_csv('./bonding_delta.csv', sep = ';').to_numpy()#importa os valores da energia de ligação da estrutura delta
atomic_mass = pd.read_csv('./atomic_mass.csv', sep = ';').to_numpy()#importa os valores da energia de ligação da estrutura delta
composition_vetor = np.zeros(len(hi_values)) #cria um vetor de zeros com a mesma quantidade de coordenadas que entradas de elementos na matriz hi_values




################# Functions

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

def Hplat(H1,c_H1,H2,c_H2,T):
    return float(H2[T][np.where(H2[T]==c_H2)[0],1] - H1[T][np.where(H1[T]==c_H1)[0],1]) / (c_H2 - c_H1)

def Splat(S1,c_H1,S2,c_H2,T):
    return float(S2[T][np.where(S2[T]==c_H2)[0],1] - S1[T][np.where(S1[T]==c_H1)[0],1]) / (c_H2 - c_H1)


def Hconvert(p,mass,cH):
    p[np.where(p==cH)[0][0],0] = cH * 1.01 *100/ (mass+ cH*1.01)
    return

#Remake the legends of PCT diagram with T only
def search (lista, valor):
    return [(lista.index(x), x.index(valor)) for x in lista if valor in x]

def dif(x,y):
    return np.absolute(y-x)




class Phase: 
        
    def __init__(self,name):
        self.name = name
        self.alloy = alloy
        self.T = []
        self.H_M = []
        self.h_m = []
        self.r = []
        self.theta = []
        self.cH_step = []
        self.cH = []
        self.wtH = []
        self.temperatures = []
        self.cH_limit = []
        self.G = {}
        self.H = {}
        self.S = {}
        self.S_c = {}
        self.mu_H = {}
        self.mu_M = {}
        self.cH_plateau_start = {}
        self.cH_plateau_end = {}
        self.pressure = {}
        self.plateau = {}
        
        
    @property
    def atomic_mass(self):
        self.composition_vetor = np.zeros(len(hi_values))
        for i in elements:
            self.composition_vetor[np.where(hi_values==i)[0]] = alloy_composition[i]
        self.alloy_atomic_mass=sum(self.composition_vetor * atomic_mass[0:,1])
        
    @property
    def H_M_calculation(self):  #determines the Enthalpy due the phase transition from the reference -- H_M diff from 0 only for delta phase
        
        if self.name == 'delta':
            alpha_vetor = np.zeros([len(bonding_alpha),len(bonding_alpha)])
            for i in elements:
                for j in elements:
                    alpha_vetor[np.where(bonding_alpha==i)[0],np.where(bonding_alpha==j)[0]]=alloy_composition[i]*alloy_composition[j]

            average_bond_alpha = alpha_vetor * bonding_alpha[0:,1:]
            average_bond_alpha = sum(sum(average_bond_alpha))

            E_total_alpha = average_bond_alpha * 4


            delta_vetor = np.zeros([len(bonding_delta),len(bonding_delta)])

            for i in elements:
                for j in elements:
                    delta_vetor[np.where(bonding_delta==i)[0],np.where(bonding_delta==j)[0]]=alloy_composition[i]*alloy_composition[j]

            average_bond_delta = delta_vetor * bonding_delta[0:,1:]
            average_bond_delta = sum(sum(average_bond_delta))
            E_total_delta = average_bond_delta * 6

            H_delta = E_total_delta - E_total_alpha
            self.H_M = H_delta
        else:
            self.H_M = 0
        
    
    
    @property
    def h_m_calculation(self):
                
        if self.name == 'alpha':
            self.h_m = sum(self.composition_vetor * hi_values[0:,1])
            
        if self.name == 'beta':
            self.h_m = sum(self.composition_vetor * hi_values[0:,2])     
            
        if self.name == 'delta':
            self.h_m = sum(self.composition_vetor * hi_values[0:,3])
    
    @property
    def set_cH_limit(self):   
        self.cH_limit = self.theta/self.r
    
             
    @property
    def set_cH(self):
        c_H = np.zeros(len(range(1,int(self.cH_limit/self.cH_step ))))
        for i in range(1,int(self.cH_limit/self.cH_step)):
            c_H[i-1] = truncate(self.cH_step*i,4)
        self.cH = c_H
        self.wtH = self.cH* 1.01 *100/(self.alloy_atomic_mass+ self.cH*1.01)
        
    @property
    def calculate_enthalpy(self):
        for T in self.T:
            self.H[T] = self.h_m * self.cH + self.H_M    
    
    def s0(self,T):
        if T !=0:
            t=T/1000
            A = 33.066178
            B = -11.363417
            C = 11.432816
            D = -2.772874
            E = -0.158558
            G = 172.707974
            S0 = A*ln(t) + B*t + C*(t**2)/2 + D*(t**3)/3 - E/(2*t**2) + G # Hydrogen standard entropy J/mol of H2- NIST
            S0 = S0/1000 #kJ/mol of H2
            return S0
        else:
            return 0

    
    @property
    def calculate_entropy(self):
        for T in self.T:
            if T!=0:
                aux1 = (self.cH)/(self.theta - ((self.r-1)*self.cH))
                aux2 = (self.theta-(self.r*self.cH))/(self.theta - ((self.r-1)*self.cH))
                self.S_c[T] = -R * (self.cH* ln(aux1) +  (self.theta- (self.r*self.cH)) * ln(aux2)) 
                self.S[T] = self.S_c[T] - ((self.cH * self.s0(T))/2)
            else:
                self.S_c[T] = np.zeros(len(self.cH))
                self.S[T] = self.S_c[T]
            
    @property
    def calculate_Gibbs(self):
        for T in self.T:
            if T!= 0:
                self.G[T] = self.H[T] - T* self.S[T]
            if T == 0:
                self.G[T] = self.H[T]
                
    @property            
    def calculate_mu_H(self):
        for T in self.T:
            aux1 = self.cH*(self.theta - (self.cH*(self.r-1)))**(self.r-1)
            aux2 = (self.theta - (self.cH*self.r))**self.r
            self.mu_H[T] = self.h_m - T * (- R* ln(aux1/aux2)- self.s0(T)/2)
            
            
    @property            
    def calculate_mu_M(self):
        for T in self.T:
            self.mu_M[T] = self.G[T] - self.cH * self.mu_H[T]




elements = []
composition = []
alloy_composition ={}
stop = False
R = 0.00831446261815324  #kJ/mol

##########Input################
print(f'The elements presents in the model are:')
print(*hi_values[0:,0].astype(str), sep= ", " )

while (stop!= True):
    y = input("Enter with an element: ")
    if str.lower(y) in np.char.lower(hi_values[0:,0].astype(str)): #verify if the element is in the model
        elements.append(y)  #add an element to the alloy
        x = float(input("Atomic fraction of element {}: " .format(y)))
        composition.append(x) #add the element to end of the composition list
        if y in alloy_composition:
            alloy_composition[y]= alloy_composition[y] + x
        else:
            alloy_composition[y]=x
    else:
        print('This element can not be used in this model yet.')
    z = input("Would you like to add another element? Y/n ")
    if z == "n":
        stop = True
        print(f'Composition of the input: {sum(composition)}')
for i in alloy_composition:
    alloy_composition[i] = alloy_composition[i]/sum(composition)  
    
    
 ########## Alloy input ##############
alloy_name = []
for element in alloy_composition.keys():
    alloy_name.append(element)
    alloy_name.append(str(truncate(alloy_composition[element],3)))

separator = ''
alloy = separator.join(alloy_name)
print(f"The normalized alloy inserted is: {alloy}")   





############# Input Temperatures ##################
Temperatures = []
stop = False
while (stop!= True):
    y = float(input("Enter a temperature in degrees Celsius for calculations: "))
    y= y + 273.15
    Temperatures.append(y) 
    z = input("Would you like to add another temperature ? Y/n ")
    if z == "n":
        stop = True
Temperatures.sort() #sort the temperatures in ascending order         


###########Input step size

ch = input("Enter the hydrogen composition step size (Default:0.00025) \ninput:")
if not ch or ch == None:
    cH_step = 0.00025
else:
    cH_step = float(ch)
    
print("Enter the theta and r entropy parameters for Alpha, Beta and Delta phases in the following order: \ntheta_alpha,r_alpha,theta_beta,r_beta,theta_delta,r_delta")
rt = (input("Default: 6,7,4,2,2,1 \ninput:"))
if not rt or rt == None:
    rt = (6,7,4,2,2,1)
if type(rt) == str:
    rt = rt.split(',') 
        




print('Your calculation has started')
begin = time.time()
#defining phases properties
############ Alpha ##############
alpha = Phase('alpha')                           #name the phase: alpha, beta or delta
alpha.atomic_mass
alpha.H_M_calculation                            #determine the enthalpy of phase transition from the reference
alpha.h_m_calculation                            #determine the enthalpy due the hyrogen on the alloy
alpha.theta = int(rt[0])                              #configurational entropy parameter: quantity of interstitial sites per metal atom
alpha.r = int(rt[1])                                  #configurational entropy parameter: site blocking effect
alpha.set_cH_limit                               #set the Hydrogen composition limit for the phase
alpha.cH_step = cH_step                          #set the step variation in the Hydrogen composition for the termodynamic calculation
alpha.set_cH                                     #set the Hydrogen compositions for the phase       
alpha.T = Temperatures                           #set the temeperatures for calculations
alpha.calculate_enthalpy                         #calculate the enthalpy of the alloy for all temperatures seted
alpha.calculate_entropy                          #calculate the entropy of the alloy for all temperatures seted
alpha.calculate_Gibbs                            #calculate the Gibbs Free Energy of the alloy for all temperatures seted
alpha.calculate_mu_H                             #calculate the Hydrogen chemical potential for all temperatures seted
alpha.calculate_mu_M                             #calculate the Metal chemical potential for all temperatures seted

########## Beta ##############
beta = Phase('beta')                            #name the phase: alpha, beta or delta
beta.atomic_mass
beta.H_M_calculation                            #determine the enthalpy of phase transition from the reference
beta.h_m_calculation                            #determine the enthalpy due the hyrogen on the alloy
beta.theta = int(rt[2])                              #configurational entropy parameter: quantity of interstitial sites per metal atom
beta.r = int(rt[3])                                  #configurational entropy parameter: site blocking effect
beta.set_cH_limit                               #set the Hydrogen composition limit for the phase
beta.cH_step = cH_step                          #set the step variation in the Hydrogen composition for the termodynamic calculation
beta.set_cH                                     #set the Hydrogen compositions for the phase       
beta.T = Temperatures                           #set the temeperatures for calculations
beta.calculate_enthalpy                         #calculate the enthalpy of the alloy for all temperatures seted
beta.calculate_entropy                          #calculate the entropy of the alloy for all temperatures seted
beta.calculate_Gibbs                            #calculate the Gibbs Free Energy of the alloy for all temperatures seted
beta.calculate_mu_H                             #calculate the Hydrogen chemical potential for all temperatures seted
beta.calculate_mu_M                             #calculate the Metal chemical potential for all temperatures seted

######### Delta ##############
delta = Phase('delta')                           #name the phase: alpha, beta or delta
delta.atomic_mass
delta.H_M_calculation                            #determine the enthalpy of phase transition from the reference
delta.h_m_calculation                            #determine the enthalpy due the hyrogen on the alloy
delta.theta = int(rt[4])                              #configurational entropy parameter: quantity of interstitial sites per metal atom
delta.r = int(rt[5])                                  #configurational entropy parameter: site blocking effect
delta.set_cH_limit                               #set the Hydrogen composition limit for the phase
delta.cH_step = cH_step                          #set the step variation in the Hydrogen composition for the termodynamic calculation
delta.set_cH                                     #set the Hydrogen compositions for the phase       
delta.T = Temperatures                           #set the temeperatures for calculations
delta.calculate_enthalpy                         #calculate the enthalpy of the alloy for all temperatures seted
delta.calculate_entropy                          #calculate the entropy of the alloy for all temperatures seted
delta.calculate_Gibbs                            #calculate the Gibbs Free Energy of the alloy for all temperatures seted
delta.calculate_mu_H                             #calculate the Hydrogen chemical potential for all temperatures seted
delta.calculate_mu_M                             #calculate the Metal chemical potential for all temperatures seted



structures = [alpha,beta,delta]     #set the structures for the termodynamic calculations




############### Search possible equilibriums compositions (condition of hydrogen and metal chemical potentials are the same in different phases)
#try:
saveab = {}
savead = {}
savebd = {}
converge = 0.01

for T in Temperatures:
    saveab[T] = np.zeros(2)
    savead[T] = np.zeros(2)
    savebd[T] = np.zeros(2)   


    for i in range(len(alpha.mu_H[T])):                              #Search the equilibrium between alpha-beta phases
        diff_array_H = np.absolute(beta.mu_H[T] - alpha.mu_H[T][i])
        crit_1 = dif(alpha.mu_H[T][i],beta.mu_H[T][diff_array_H.argmin()])
        crit_2 = dif(alpha.mu_M[T][i],beta.mu_M[T][diff_array_H.argmin()])

        if (min(diff_array_H) < converge)  and (crit_1<converge) and (crit_2<converge) and (alpha.cH[i] < beta.cH[diff_array_H.argmin()]):
            saveab[T][0]=alpha.cH[i]
            saveab[T][1]=beta.cH[diff_array_H.argmin()]
            break

        diff_array_M = np.absolute(beta.mu_M[T] - alpha.mu_M[T][i])
        crit_3 = dif(alpha.mu_H[T][i],beta.mu_H[T][diff_array_M.argmin()])
        crit_4 = dif(alpha.mu_M[T][i],beta.mu_M[T][diff_array_M.argmin()])

        if (min(diff_array_M) < converge)  and (crit_3<converge) and (crit_4<converge) and (alpha.cH[i] < beta.cH[diff_array_M.argmin()]):
            saveab[T][0]= alpha.cH[i]
            saveab[T][1]= beta.cH[diff_array_M.argmin()]
            if any(saveab[T]):
                break
                
                
    for i in range(len(beta.mu_H[T])):                               #Search the equilibrium between beta-delta phases
        diff_array_H = np.absolute(delta.mu_H[T] - beta.mu_H[T][i])
        crit_1 = dif(beta.mu_H[T][i],delta.mu_H[T][diff_array_H.argmin()])
        crit_2 = dif(beta.mu_M[T][i],delta.mu_M[T][diff_array_H.argmin()])

        if (min(diff_array_H) < converge)  and (crit_1<converge) and (crit_2<converge) and (beta.cH[i] < delta.cH[diff_array_H.argmin()]):
            savebd[T][0]= beta.cH[i]
            savebd[T][1]= delta.cH[diff_array_H.argmin()]
            break

        diff_array_M = np.absolute(delta.mu_M[T] - beta.mu_M[T][i])
        crit_3 = dif(beta.mu_H[T][i],delta.mu_H[T][diff_array_M.argmin()])
        crit_4 = dif(beta.mu_M[T][i],delta.mu_M[T][diff_array_M.argmin()])

        if (min(diff_array_M) < converge)  and (crit_3<converge) and (crit_4<converge) and (beta.cH[i] < delta.cH[diff_array_M.argmin()]):
            savebd[T][0]= beta.cH[i]
            savebd[T][1]= delta.cH[diff_array_M.argmin()]
            break  

    for i in range(len(alpha.mu_H[T])):                               #Search the equilibrium between alpha-delta phases
        diff_array_H = np.absolute(delta.mu_H[T] - alpha.mu_H[T][i])
        crit_1 = dif(alpha.mu_H[T][i],delta.mu_H[T][diff_array_H.argmin()])
        crit_2 = dif(alpha.mu_M[T][i],delta.mu_M[T][diff_array_H.argmin()])

        if (min(diff_array_H) < converge)  and (crit_1<converge) and (crit_2<converge) and (alpha.cH[i] < delta.cH[diff_array_H.argmin()]):
            savead[T][0]= alpha.cH[i]
            savead[T][1]= delta.cH[diff_array_H.argmin()]
            break           

        diff_array_M = np.absolute(delta.mu_M[T] - alpha.mu_M[T][i])
        crit_3 = dif(alpha.mu_H[T][i],delta.mu_H[T][diff_array_M.argmin()])
        crit_4 = dif(alpha.mu_M[T][i],delta.mu_M[T][diff_array_M.argmin()])

        if (min(diff_array_M) < converge)  and (crit_3<converge) and (crit_4<converge) and (alpha.cH[i] < delta.cH[diff_array_M.argmin()]):
            savead[T][0]= alpha.cH[i]
            savead[T][1]= delta.cH[diff_array_M.argmin()]
            break  





equilibrium_real = {}    
equilibrium_all = {}  
verify_equilibrium = {}
for T in Temperatures:
    if T != 0:
        equilibrium_all[T] = {'ab': saveab[T], 'ad' : savead[T], 'bd' : savebd[T]}
        
        if delta.H_M >0:
            if any(saveab[T]):
                if any(savead[T]) and (saveab[T][0]>savead[T][0]):
                    equilibrium_real[T] = 'ad'
                    alpha.cH_plateau_start[T] = savead[T][0]
                    delta.cH_plateau_end[T] = savead[T][1]
                    beta.cH_plateau_end[T] = None
                    beta.cH_plateau_start[T] = None 
                if not any(savead[T]) or (saveab[T][0]<savead[T][0]):
                    equilibrium_real[T] = 'ab'
                    alpha.cH_plateau_start[T] = saveab[T][0]
                    beta.cH_plateau_end[T] = saveab[T][1]

                    if any(savebd[T]):
                        beta.cH_plateau_start[T] = savebd[T][0]
                        delta.cH_plateau_end[T] = savebd[T][1]
                    if not any(savebd[T]):
                        search_cH_theta_r = dif(beta.theta/beta.r,beta.cH)
                        beta.cH_plateau_start[T] = beta.cH[np.where(search_cH_theta_r == min(search_cH_theta_r))]
                        delta.cH_plateau_end[T] = delta.cH[-1]
                        if beta.cH_plateau_start[T] == delta.cH_plateau_end[T]:
                            print('You must change the values of Theta and r for Beta phase')
                            break
                            
            if any(savead[T]) and not any(saveab[T]):
                equilibrium_real[T] = 'ad'
                alpha.cH_plateau_start[T] = savead[T][0]
                delta.cH_plateau_end[T] = savead[T][1]
                beta.cH_plateau_end[T] = None
                beta.cH_plateau_start[T] = None     
       
            if not (any(saveab[T]) or any(savead[T])):
                print('Equilibrium not found, you must check the input parameters')
                break
                
                
        if delta.H_M <0:
            if any(saveab[T]):
                equilibrium_real[T] = 'ab'
                alpha.cH_plateau_start[T] = saveab[T][0]
                beta.cH_plateau_end[T] = saveab[T][1]
                beta.cH_plateau_start[T] = None
            if not any(saveab[T]):
                equilibrium_real[T] = 'Non-equilibrium'
                alpha.cH_plateau_start[T] = None
                beta.cH_plateau_end[T] = None
                beta.cH_plateau_start[T] = None
            
            




plat_ab ={}
plat_bd = {}
plat_ad = {}
H_plat_ab = {}
S_plat_ab = {}
H_plat_bd = {}
S_plat_bd = {}
H_plat_ad = {}
S_plat_ad = {}
P_0 = 1

for T in Temperatures:
    if T!= 0:

        if delta.H_M <0:    #condition to delta-phase is more stable than alpha or beta
            if not alpha.cH_plateau_start:
                alpha.pressure[T] = np.zeros(len(alpha.cH))
                for c_h_index in len(alpha.cH):
                    mu_alpha = float(alpha.mu_H[c_h_index])
                    lnpH = 2* mu_alpha/(R*T)
                    peq = P_0*exp(lnpH)
                    alpha.pressure[T] = peq 
                    
            if beta.cH_plateau_start[T] == None:
                p_alpha_temp = np.zeros(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])
                for c_h_index in range(int(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])):
                    mu_alpha = float(alpha.mu_H[T][c_h_index])
                    lnpH = 2* mu_alpha/(R*T)
                    peq = P_0*exp(lnpH)
                    p_alpha_temp[c_h_index] = peq 
                alpha.pressure[T] = p_alpha_temp

                ar_size = len(beta.cH)-np.where(beta.cH == beta.cH_plateau_end[T])[0]
                p_beta_temp = np.zeros(ar_size)                
                for c_h_index in range(int(np.where(beta.cH == beta.cH_plateau_end[T])[0]),len(beta.cH)):
                    mu_beta = float(beta.mu_H[T][c_h_index])
                    lnpH = 2* mu_beta/(R*T)
                    peq = P_0*exp(lnpH)
                    p_beta_temp[c_h_index - int(np.where(beta.cH == beta.cH_plateau_end[T])[0])] = peq
                beta.pressure[T] = p_beta_temp
                
                ############### Determine plateau alpha-beta phases

                H_plat_ab[T] =  (beta.H[T][np.where(beta.cH == beta.cH_plateau_end[T])] - alpha.H[T][np.where(alpha.cH == alpha.cH_plateau_start[T])])/ (beta.cH_plateau_end[T] - alpha.cH_plateau_start[T])  
                S_plat_ab[T] =  (beta.S[T][np.where(beta.cH == beta.cH_plateau_end[T])] - alpha.S[T][np.where(alpha.cH == alpha.cH_plateau_start[T])])/ (beta.cH_plateau_end[T] - alpha.cH_plateau_start[T])
                plateau_ab = exp(2* ((H_plat_ab[T]/(R*T) - (S_plat_ab[T]/R))))
                plat_ab[T] = np.full(2,plateau_ab)
                alpha.plateau[T] = plat_ab[T][0]
                
                
        if delta.H_M >0:
            if alpha.cH_plateau_start[T] != None:

                p_alpha_temp = np.zeros(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])        
                for c_h_index in range(int(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])):
                    mu_alpha = float(alpha.mu_H[T][c_h_index])
                    lnpH = 2* mu_alpha/(R*T)
                    peq = P_0*exp(lnpH)
                    p_alpha_temp[c_h_index] = peq 
                alpha.pressure[T] = p_alpha_temp


            if beta.cH_plateau_start[T] != None:


                ############### Determine plateau alpha-beta phases

                H_plat_ab[T] =  (beta.H[T][np.where(beta.cH == beta.cH_plateau_end[T])] - alpha.H[T][np.where(alpha.cH == alpha.cH_plateau_start[T])])/ (beta.cH_plateau_end[T] - alpha.cH_plateau_start[T])  
                S_plat_ab[T] =  (beta.S[T][np.where(beta.cH == beta.cH_plateau_end[T])] - alpha.S[T][np.where(alpha.cH == alpha.cH_plateau_start[T])])/ (beta.cH_plateau_end[T] - alpha.cH_plateau_start[T])
                plateau_ab = exp(2* ((H_plat_ab[T]/(R*T) - (S_plat_ab[T]/R))))
                plat_ab[T] = np.full(2,plateau_ab)
                alpha.plateau[T] = plat_ab[T][0]
                
                ########## Determine Plateau beta - delta phases
                H_plat_bd[T] =  (delta.H[T][np.where(delta.cH == delta.cH_plateau_end[T])] - beta.H[T][np.where(beta.cH == beta.cH_plateau_start[T])])/ (delta.cH_plateau_end[T] - beta.cH_plateau_start[T])  
                S_plat_bd[T] =  (delta.S[T][np.where(delta.cH == delta.cH_plateau_end[T])] - beta.S[T][np.where(beta.cH == beta.cH_plateau_start[T])])/ (delta.cH_plateau_end[T] - beta.cH_plateau_start[T])
                plateau_bd = exp(2* ((H_plat_bd[T]/(R*T) - (S_plat_bd[T]/R))))
                plat_bd[T] = np.full(2,plateau_bd)
                beta.plateau[T]= plat_bd[T][0]

                #############Determine pressure of beta phase
                ar_size = np.where(beta.cH == beta.cH_plateau_start[T])[0]-np.where(beta.cH == beta.cH_plateau_end[T])[0] #+1
                p_beta_temp = np.zeros(ar_size)
                P_0 = 1        
                for c_h_index in range(int(np.where(beta.cH == beta.cH_plateau_end[T])[0]),int(np.where(beta.cH == beta.cH_plateau_start[T])[0])):#+1):
                    mu_beta = float(beta.mu_H[T][c_h_index])
                    lnpH = 2* mu_beta/(R*T)
                    peq = P_0*exp(lnpH)
                    if peq < beta.plateau[T]: #When we assume that on the beta phase H composition theta/r the saturation occurs, for cH close to theta/r the peq explodes to infinity
                        p_beta_temp[c_h_index - int(np.where(beta.cH == beta.cH_plateau_end[T])[0])] = peq
                    else:
                        p_beta_temp[c_h_index - int(np.where(beta.cH == beta.cH_plateau_end[T])[0])] = beta.plateau[T]
                beta.pressure[T] = p_beta_temp


            if beta.cH_plateau_start[T] == None:
                beta.pressure[T] = None
                ########## Determine Plateau alpha - delta phases
                H_plat_ad[T] =  (delta.H[T][np.where(delta.cH == delta.cH_plateau_end[T])] - alpha.H[T][np.where(alpha.cH == alpha.cH_plateau_start[T])])/ (delta.cH_plateau_end[T] - alpha.cH_plateau_start[T])  
                S_plat_ad[T] =  (delta.S[T][np.where(delta.cH == delta.cH_plateau_end[T])] - alpha.S[T][np.where(alpha.cH == alpha.cH_plateau_start[T])])/ (delta.cH_plateau_end[T] - alpha.cH_plateau_start[T])
                plateau_ad = exp(2* ((H_plat_ad[T]/(R*T) - (S_plat_ad[T]/R))))
                plat_ad[T] = np.full(2,plateau_ad)
                delta.plateau[T] = plat_ad[T][0]

            if delta.cH_plateau_end[T] != None:
                ar_size = len(delta.cH) - np.where(delta.cH == delta.cH_plateau_end[T])[0]
                p_delta_temp = np.zeros(ar_size)
                P_0 = 1        
                for c_h_index in range(int(np.where(delta.cH == delta.cH_plateau_end[T])[0]),len(delta.cH)):
                    mu_delta = float(delta.mu_H[T][c_h_index])
                    lnpH = 2* mu_delta/(R*T)
                    peq = P_0*exp(lnpH)
                    p_delta_temp[c_h_index - int(np.where(delta.cH == delta.cH_plateau_end[T])[0])] = peq
                delta.pressure[T] = p_delta_temp







#########################################
#  Saving and ploting files and data   #  
######################################### 
if not os.path.isdir(f'./{alloy}_RM'): #Verify if the folder exists
    os.mkdir(f'./{alloy}_RM')          #creates the folder is it not exists
    
dnow = datetime.now()        
my_path = os.path.abspath(f'./{alloy}_RM')
my_file = f'./{alloy}.txt'
if not os.path.isfile(os.path.join(my_path, my_file)):
    file = open(os.path.join(my_path, my_file), 'w+')
    file.close()
file = open(os.path.join(my_path, my_file),'a')
file.write(f"Alloy: {alloy} | {dnow.strftime('%d/%m/%Y %H:%M')} | \n")
file.write(f"hM[alpha]= {truncate(alpha.h_m,3)} | hM[beta]= {truncate(beta.h_m,3)} | hM[delta]= {truncate(delta.h_m,3)} |\n")
file.write(f"HM[alpha] = {truncate(alpha.H_M,3)}| HM[beta] = {truncate(beta.H_M,3)} | HM[delta] = {truncate(delta.H_M,3)} \n")
file.write("Values for Theta and r is in the format: structure(theta,r)\n")
file.write(f"alpha:{alpha.theta,alpha.r}, beta:{beta.theta,beta.r}, delta:{delta.theta,delta.r}\n")
file.write(f"cH step - alpha:{alpha.cH_step}, beta: {beta.cH_step}, delta: {delta.cH_step}\n")
file.write(f"------------------------------------------------------- \n")
file.write(f"####################################################### \n")
file.write(f"------------------------------------------------------- \n")
for T in Temperatures:
    file.write(f"Temperature = {T-273.15} \u00B0C | {T} K |\n")
    if equilibrium_real[T] == 'ab':
        file.write(f"Plateau alfa-beta \n")
        file.write(f"c_H[alpha] = {alpha.cH_plateau_start[T]} \n")
        file.write(f"c_H[beta] = {beta.cH_plateau_end[T]} \n")
        file.write(f"%wt H[alpha] = {float(alpha.wtH[np.where(alpha.cH ==alpha.cH_plateau_start[T])])} \n")
        file.write(f"%wt H[beta] = {float(beta.wtH[np.where(beta.cH ==beta.cH_plateau_end[T])])} \n")
        file.write(f"Pressure plateau = {alpha.plateau[T]} [atm] \n")
        file.write(f"Enthalpy plateau = {float(H_plat_ab[T])} kJ/mol H \n")
        file.write(f"Entropy plateau = {float(S_plat_ab[T])} kJ/mol H \n")
        file.write(f"____________________________________________________ \n")
        if beta.cH_plateau_start[T] != None:
            file.write(f"Plateau beta-delta \n")
            file.write(f"c_H[beta] = {beta.cH_plateau_start[T]} \n")
            file.write(f"c_H[delta] = {delta.cH_plateau_end[T]} \n")
            file.write(f"%wt H[beta] = {float(beta.wtH[np.where(beta.cH ==beta.cH_plateau_start[T])])} \n")
            file.write(f"%wt H[delta] = {float(delta.wtH[np.where(delta.cH ==delta.cH_plateau_end[T])])} \n")
            file.write(f"Pressure plateau = {beta.plateau[T]} [atm] \n")
            file.write(f"Enthalpy plateau = {float(H_plat_bd[T])} kJ/mol H \n")
            file.write(f"Entropy plateau = {float(S_plat_bd[T])} kJ/mol H \n")

    if equilibrium_real[T] == 'ad':
        file.write(f"Plateau alfa-delta \n")
        file.write(f"c_H[alpha] = {alpha.cH_plateau_start[T]} \n")
        file.write(f"c_H[delta] = {delta.cH_plateau_end[T]} \n")
        file.write(f"%wt H[alpha] = {float(alpha.wtH[np.where(alpha.cH ==alpha.cH_plateau_start[T])])} \n")
        file.write(f"%wt H[delta] = {float(delta.wtH[np.where(delta.cH ==delta.cH_plateau_end[T])])} \n")
        file.write(f"Pressure plateau = {delta.plateau[T]} [atm] \n")
        file.write(f"Enthalpy plateau = {float(H_plat_ad[T])} kJ/mol H \n")
        file.write(f"Entropy plateau = {float(S_plat_ad[T])} kJ/mol H \n")        

    file.write(f"------------------------------------------------------- \n")
    file.write(f"####################################################### \n")
    file.write(f"------------------------------------------------------- \n")
file.close()
########################################################################################




properties = ['H','S','S_c','G']
properties2 = ['mu_H','mu_M']

for T in Temperatures:
    if not os.path.isdir(f'./{alloy}_RM/{T - 273.15}'): 
        os.mkdir(f'./{alloy}_RM/{T - 273.15}')
    my_path = os.path.abspath(alloy+ f"_RM/{T - 273.15}") # Figures out the absolute path for you in case your working directory moves around.

    for prop in properties:
     
    
        my_file = f'{alloy}_{prop}_T_{T-273.15}_RM.png'
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(alpha.cH, alpha.__dict__[prop][T], 'b', label = f'$\Delta {prop}$' + r'$^{\alpha}$')
        axes.plot(beta.cH, beta.__dict__[prop][T], 'orange', label = f'$\Delta {prop}$' + r'$^{\beta}$')
        axes.plot(delta.cH, delta.__dict__[prop][T], 'g', label = f'$\Delta {prop}$' + r'$^{\delta}$')
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(f"$\Delta {prop}$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')

        ####################################### Plot %wt H ###############################        
        my_file = f'{alloy}_{prop}_T_{T-273.15}_wtH_RM.png'
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(alpha.wtH, alpha.__dict__[prop][T], 'b', label = f'$\Delta {prop}$' + r'$^{\alpha}$')
        axes.plot(beta.wtH, beta.__dict__[prop][T], 'orange', label = f'$\Delta {prop}$' + r'$^{\beta}$')
        axes.plot(delta.wtH, delta.__dict__[prop][T], 'g', label = f'$\Delta {prop}$' + r'$^{\delta}$')
        axes.set_xlabel('%wt H') # Notice the use of set_ to begin methods
        axes.set_ylabel(f"$\Delta {prop}$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')       
        
     
      
    for prop in properties2:
     
    
        my_file = f'{alloy}_{prop}_T_{T-273.15}_RM.png'
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(alpha.cH, alpha.__dict__[prop][T], 'b', label = f'$\Delta \{prop}$' + r'$^{\alpha}$')
        axes.plot(beta.cH, beta.__dict__[prop][T], 'orange', label = f'$\Delta \{prop}$' + r'$^{\beta}$')
        axes.plot(delta.cH, delta.__dict__[prop][T], 'g', label = f'$\Delta \{prop}$' + r'$^{\delta}$')
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(f"$\Delta \{prop}$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')

        ####################################### Plot %wt H ###############################          
        my_file = f'{alloy}_{prop}_T_{T-273.15}_wtH_RM.png'
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(alpha.wtH, alpha.__dict__[prop][T], 'b', label = f'$\Delta \{prop}$' + r'$^{\alpha}$')
        axes.plot(beta.wtH, beta.__dict__[prop][T], 'orange', label = f'$\Delta \{prop}$' + r'$^{\beta}$')
        axes.plot(delta.wtH, delta.__dict__[prop][T], 'g', label = f'$\Delta \{prop}$' + r'$^{\delta}$')
        axes.set_xlabel('%wt H') # Notice the use of set_ to begin methods
        axes.set_ylabel(f"$\Delta \{prop}$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')




for T in Temperatures:
        if T!=0:
            my_path = os.path.abspath(alloy+ f"_RM/{T - 273.15}")
            my_file = f'{alloy}_PCI_T_{T-273.15}_RM.png'
            fig = plt.figure()

            try:
                sup_range_alpha = int(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])
            except:
                print(f'The alloy does not shows other phases at {T- 273.15} \u00B0C')
            try:
                low_range_beta = int(np.where(beta.cH == beta.cH_plateau_end[T])[0])
                if delta.H_M >0:
                    sup_range_beta = int(np.where(beta.cH == beta.cH_plateau_start[T])[0])
            except:
                print(f'The alloy does not shows beta-phase at {T- 273.15} \u00B0C')
                
            if delta.H_M >0:
                try:
                    low_range_delta = int(np.where(delta.cH == delta.cH_plateau_end[T])[0])
                    sup_range_delta = len(delta.cH)
                except:
                    print(f'The alloy does not shows delta-phase at {T- 273.15} \u00B0C')
            # Add set of axes to figure
            axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

            # Plot on that set of axes
            axes.plot(alpha.cH[0:sup_range_alpha], alpha.pressure[T][0:sup_range_alpha], 'b', label = r'$\alpha$')
            
            if delta.H_M <0: 
                axes.plot([alpha.cH[sup_range_alpha],beta.cH[low_range_beta]], plat_ab[T], 'red')
                axes.plot(beta.cH[low_range_beta:], beta.pressure[T], 'orange', label = r'$\beta$')

            if beta.cH_plateau_start[T] != None:
                axes.plot([alpha.cH[sup_range_alpha],beta.cH[low_range_beta]], plat_ab[T], 'red')
                axes.plot(beta.cH[low_range_beta:sup_range_beta], beta.pressure[T], 'orange', label = r'$\beta$')
                axes.plot([beta.cH[sup_range_beta-1],delta.cH[low_range_delta]], plat_bd[T], 'red')


            if (beta.cH_plateau_start[T] == None) and (delta.H_M > 0):
                axes.plot([alpha.cH[sup_range_alpha],delta.cH[low_range_delta]], plat_ad[T], 'red')

            if delta.H_M > 0:    
                axes.plot(delta.cH[low_range_delta:sup_range_delta], delta.pressure[T], 'g', label = r'$\delta$')

            plt.yscale("log")
            axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
            axes.set_ylabel(r"$Pressure$ [atm]")
            axes.set_title(f'T = {T- 273.15} \u00B0C')
            axes.legend(loc = 0)
            fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
            
################################################ Plot %wt H ###############################
            my_path = os.path.abspath(alloy+ f"_RM/{T - 273.15}")
            my_file = f'{alloy}_PCI_T_{T-273.15}_RM.png'
            fig = plt.figure()

            try:
                sup_range_alpha = int(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])
            except:
                print(f'The alloy does not shows other phases at {T- 273.15} \u00B0C')
            try:
                low_range_beta = int(np.where(beta.cH == beta.cH_plateau_end[T])[0])
                if delta.H_M >0:
                    sup_range_beta = int(np.where(beta.cH == beta.cH_plateau_start[T])[0])
            except:
                print(f'The alloy does not shows beta-phase at {T- 273.15} \u00B0C')
                
            if delta.H_M >0:
                try:
                    low_range_delta = int(np.where(delta.cH == delta.cH_plateau_end[T])[0])
                    sup_range_delta = len(delta.cH)
                except:
                    print(f'The alloy does not shows delta-phase at {T- 273.15} \u00B0C')
            # Add set of axes to figure
            axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

            # Plot on that set of axes
            axes.plot(alpha.wtH[0:sup_range_alpha], alpha.pressure[T][0:sup_range_alpha], 'b', label = r'$\alpha$')
            
            if delta.H_M <0: 
                axes.plot([alpha.wtH[sup_range_alpha],beta.wtH[low_range_beta]], plat_ab[T], 'red')
                axes.plot(beta.wtH[low_range_beta:], beta.pressure[T], 'orange', label = r'$\beta$')

            if beta.cH_plateau_start[T] != None:
                axes.plot([alpha.wtH[sup_range_alpha],beta.wtH[low_range_beta]], plat_ab[T], 'red')
                axes.plot(beta.wtH[low_range_beta:sup_range_beta], beta.pressure[T], 'orange', label = r'$\beta$')
                axes.plot([beta.wtH[sup_range_beta-1],delta.wtH[low_range_delta]], plat_bd[T], 'red')


            if (beta.cH_plateau_start[T] == None) and (delta.H_M > 0):
                axes.plot([alpha.wtH[sup_range_alpha],delta.wtH[low_range_delta]], plat_ad[T], 'red')

            if delta.H_M > 0:    
                axes.plot(delta.wtH[low_range_delta:sup_range_delta], delta.pressure[T], 'g', label = r'$\delta$')

            plt.yscale("log")
            axes.set_xlabel('%wt H') # Notice the use of set_ to begin methods
            axes.set_ylabel(r"$Pressure$ [atm]")
            axes.set_title(f'T = {T- 273.15} \u00B0C')
            axes.legend(loc = 0)
            fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')




##################### PCT plot ##################################3
  
my_path = os.path.abspath(alloy+ f"_RM")
my_file = f'{alloy}_PCT_RM.png'
fig = plt.figure()
# Add set of axes to figure
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)


cmap = plt.get_cmap("Set1")
i = 0
for T in Temperatures: 
    if T!=0: 
                
        try:
            sup_range_alpha = int(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])
        except:
            print(f'The alloy does not shows other phases at {T- 273.15} \u00B0C')
        try:
            low_range_beta = int(np.where(beta.cH == beta.cH_plateau_end[T])[0])
            if delta.H_M >0:
                sup_range_beta = int(np.where(beta.cH == beta.cH_plateau_start[T])[0])
        except:
            print(f'The alloy does not shows beta-phase at {T- 273.15} \u00B0C')

        if delta.H_M >0:
            try:
                low_range_delta = int(np.where(delta.cH == delta.cH_plateau_end[T])[0])
                sup_range_delta = len(delta.cH)
            except:
                print(f'The alloy does not shows delta-phase at {T- 273.15} \u00B0C')                
                
                
                
                
        # Add set of axes to figure
        

        # Plot on that set of axes
        axes.plot(alpha.cH[0:sup_range_alpha], alpha.pressure[T][0:sup_range_alpha], color = cmap(i), label = f'T = {T-273.15} \u00B0C')

        if delta.H_M <0: 
            axes.plot([alpha.cH[sup_range_alpha],beta.cH[low_range_beta]], plat_ab[T], color = cmap(i))
            axes.plot(beta.cH[low_range_beta:], beta.pressure[T], color = cmap(i))

        if beta.cH_plateau_start[T] != None:
            axes.plot([alpha.cH[sup_range_alpha],beta.cH[low_range_beta]], plat_ab[T], color = cmap(i))
            axes.plot(beta.cH[low_range_beta:sup_range_beta], beta.pressure[T], color = cmap(i))
            axes.plot([beta.cH[sup_range_beta],delta.cH[low_range_delta]], plat_bd[T], color = cmap(i))


        if (beta.cH_plateau_start[T] == None) and (delta.H_M > 0):
            axes.plot([alpha.cH[sup_range_alpha],delta.cH[low_range_delta]], plat_ad[T], color = cmap(i))

        if delta.H_M > 0:    
            axes.plot(delta.cH[low_range_delta:sup_range_delta], delta.pressure[T], color = cmap(i))
        i += 1
plt.yscale("log")
axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
axes.set_ylabel(r"$Pressure$ [atm]")
axes.set_title(f'PCT {alloy}')
axes.legend(loc = 0, ncol=2)
fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')        


my_path = os.path.abspath(alloy+ f"_RM")
my_file = f'{alloy}_PCT_wtH_RM.png'
fig = plt.figure()
# Add set of axes to figure
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)


cmap = plt.get_cmap("Set1")
i = 0
for T in Temperatures: 
    if T!=0: 
        try:
            sup_range_alpha = int(np.where(alpha.cH == alpha.cH_plateau_start[T])[0])#+1
        except:
            print(f'The alloy does not shows other phases at {T- 273.15} \u00B0C')
        try:
            low_range_beta = int(np.where(beta.cH == beta.cH_plateau_end[T])[0])
            if delta.H_M >0:
                sup_range_beta = int(np.where(beta.cH == beta.cH_plateau_start[T])[0])#+1
        except:
            print(f'The alloy does not shows beta-phase at {T- 273.15} \u00B0C')

        if delta.H_M >0:
            try:
                low_range_delta = int(np.where(delta.cH == delta.cH_plateau_end[T])[0])
                sup_range_delta = len(delta.cH)
            except:
                print(f'The alloy does not shows delta-phase at {T- 273.15} \u00B0C')
        # Add set of axes to figure
        

        # Plot on that set of axes
        axes.plot(alpha.wtH[0:sup_range_alpha], alpha.pressure[T][0:sup_range_alpha], color = cmap(i), label = f'T = {T-273.15} \u00B0C')

        if delta.H_M <0: 
            axes.plot([alpha.wtH[sup_range_alpha],beta.wtH[low_range_beta]], plat_ab[T], color = cmap(i))
            axes.plot(beta.wtH[low_range_beta:], beta.pressure[T], color = cmap(i))

        if beta.cH_plateau_start[T] != None:
            axes.plot([alpha.wtH[sup_range_alpha],beta.wtH[low_range_beta]], plat_ab[T], color = cmap(i))
            axes.plot(beta.wtH[low_range_beta:sup_range_beta], beta.pressure[T], color = cmap(i))
            axes.plot([beta.wtH[sup_range_beta],delta.wtH[low_range_delta]], plat_bd[T], color = cmap(i))


        if (beta.cH_plateau_start[T] == None) and (delta.H_M > 0):
            axes.plot([alpha.wtH[sup_range_alpha],delta.wtH[low_range_delta]], plat_ad[T], color = cmap(i))

        if delta.H_M > 0:    
            axes.plot(delta.wtH[low_range_delta:sup_range_delta], delta.pressure[T], color = cmap(i))

        i += 1
plt.yscale("log")
axes.set_xlabel('%wt H') # Notice the use of set_ to begin methods
axes.set_ylabel(r"$Pressure$ [atm]")
axes.set_title(f'PCT {alloy}')
axes.legend(loc = 0, ncol=2)
fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')   



if len(Temperatures)>1:
    my_path = os.path.abspath(alloy+ f"_RM")
    my_file = f'{alloy}_VantHoff.png'
    fig = plt.figure()
    # Add set of axes to figure
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

    if len(beta.plateau) != False:
        axes.plot(1/np.array(list(alpha.plateau.keys())),ln(list(alpha.plateau.values())), label = r'$\alpha - \beta$ plateau')
        axes.plot(1/np.array(list(beta.plateau.keys())),ln(list(beta.plateau.values())), label = r'$\beta - \delta$ plateau')

    if len(delta.plateau) != False:
        axes.plot(1/np.array(list(delta.plateau.keys())),ln(list(delta.plateau.values())), label = r"$\alpha - \delta$ plateau")

    axes.set_ylabel('ln P$_{plateau}$')
    axes.set_xlabel(r'1/T [$K^{-1}$]')
    axes.set_title(f"{alloy} - Van't Hoff plot")
    axes.legend(loc=0)
    fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')  



properties = ['H','S','S_c','G','mu_H','mu_M']
for T in Temperatures:
    
    if not os.path.isdir(f'./{alloy}_RM/{T - 273.15}/data'):
            os.mkdir(f'./{alloy}_RM/{T - 273.15}/data')
    my_path = os.path.abspath(alloy+ f"_RM/{T - 273.15}/data") 

    
    

    
    

    for structure in structures:
        transfer =np.zeros((len(structure.cH),2))
        transfer[:,0] = structure.cH
        
        for i in properties:
        
            transfer[:,1] = structure.__dict__[i][T]  #Saving each propertie for each phase
            my_file = f'{alloy}_{i}_{structure.name}_T{T - 273.15}C_RM.txt' #structure.name[0]
            dataframe = pd.DataFrame(transfer)
            dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : f"{i}_{structure.name}"})
            dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', f"{i}_{structure.name}"])

            
############ Save wt H data ##################            
        transfer[:,0] = structure.wtH
        
        for i in properties:
        
            transfer[:,1] = structure.__dict__[i][T]  #Saving each propertie for each phase
            my_file = f'{alloy}_{i}_{structure.name}_T{T - 273.15}C_RM.txt' #structure.name[0]
            dataframe = pd.DataFrame(transfer)
            dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : f"{i}_{structure.name}"})
            dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', f"{i}_{structure.name}"])
        
       
            
end = time.time()
tempo = end - begin
print(f"The calculation takes {truncate(tempo/60,2)} minutes")     
            
 

