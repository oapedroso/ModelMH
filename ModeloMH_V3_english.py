#!/usr/bin/env python
# coding: utf-8

# In[11]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import log as ln
from numpy import exp
import os
import math
from itertools import zip_longest
from datetime import datetime  
import time

hi_values = pd.read_csv('./hi.csv', sep = ';').to_numpy() #importa os valores de hi do arquivo csv e transforma em uma matriz 
bonding_alpha = pd.read_csv('./bonding_alpha.csv', sep = ';').to_numpy() #importa os valores da energia de ligação da estrutura alpha
bonding_delta = pd.read_csv('./bonding_delta.csv', sep = ';').to_numpy()#importa os valores da energia de ligação da estrutura delta
atomic_mass = pd.read_csv('./atomic_mass.csv', sep = ';').to_numpy()#importa os valores da energia de ligação da estrutura delta
composition_vetor = np.zeros(len(hi_values)) #cria um vetor de zeros com a mesma quantidade de coordenadas que entradas de elementos na matriz hi_values

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



##########Definição de variaveis######################
elements = []
composition = []
alloy_composition ={}
stop = False
R = 0.00831446261815324

##########Entrada de dados pelo usuário################
while (stop!= True):
    y = input("Enter with an element: ")
    elements.append(y) #adiciona a entrada y no fim da lista de elementos
    x = float(input("Atomic fraction of element {}: " .format(y)))
    composition.append(x) #adiciona a entrada w no fim da lista de composições
    if y in alloy_composition:
        alloy_composition[y]= alloy_composition[y] + x
    else:
        alloy_composition[y]=x
    #if sum(composition)>1:
     #   print("A soma composição deve ser igual a 1, reinicie o programa")
    z = input("Would you like to add another element? Y/n ")
    if z == "n":
        stop = True
        print(sum(composition))
for i in alloy_composition:
    alloy_composition[i] = alloy_composition[i]/sum(composition)    


# In[12]:


##########Informa ao usuário a liga que foi inserida##############
alloy_name = []
for element in alloy_composition.keys():
    alloy_name.append(element)
    alloy_name.append(str(truncate(alloy_composition[element],3)))

separator = ''
alloy = separator.join(alloy_name)
print(f"The alloy inserted is: {alloy}")
#print(f"Com concentração:{alloy_composition}")


##########Cálculo de hi para cada estrutura ####################

for i in elements:
    composition_vetor[np.where(hi_values==i)[0]] = alloy_composition[i]

hm_alpha = composition_vetor * hi_values[0:,1]
hm_beta = composition_vetor * hi_values[0:,2]
hm_delta = composition_vetor * hi_values[0:,3]
hm_alpha= sum(hm_alpha)
hm_beta = sum(hm_beta)
hm_delta = sum(hm_delta)


#print(hm_alpha)
#print(hm_beta)
#print(hm_delta)


# In[13]:



##########Cálculo de H ############

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
H_alpha = 0
H_beta = 0
#print(H_delta)


# In[14]:



########## Input de temperatura
temperature = []
stop = False
while (stop!= True):
    y = float(input("Enter a temperature in degrees Celsius for calculations: "))
    y= y + 273.15
    temperature.append(y) #adiciona a entrada y no fim da lista de elementos
    z = input("Would you like to add another temperature ? Y/n ")
    if z == "n":
        stop = True

temperature.sort() #organiza as temperaturas em ordem crescente
##########Cálculos termodinâmicos########### 
val_S_alpha = {}
val_S_beta = {}
val_S_delta = {}
val_H_alpha = {}
val_H_beta = {}
val_H_delta = {}
val_G_alpha = {}
val_G_beta = {}
val_G_delta = {}
val_mu_M_alpha = {}
val_mu_M_beta = {}
val_mu_M_delta = {}
val_mu_H_alpha = {}
val_mu_H_beta = {}
val_mu_H_delta = {}
valoresalpha = {}
valoresbeta = {}
valoresdelta = {}
p_total = {}
p_total_wt = {}
p_plat = {}
H_plat = {}
S_plat = {}
verify_equilibrium = []
C_H = []
p_equilibrium_ad = []
p_equilibrium_ab = []
p_equilibrium_bd = []
p_notequilibrium = []
ch = input("Enter the hydrogen composition step size (Default:0.0005) \ninput:")
if not ch or ch == None:
    cH_step = 0.0005
else:
    cH_step = float(ch)
    
print("Enter the theta and r entropy parameters for Alpha, Beta and Delta phases in the following order: \ntheta_alpha,r_alpha,theta_beta,r_beta,theta_delta,r_delta")
rt = (input("Default: 6,7,4,2,2,1 \ninput:"))
if not rt or rt == None:
    rt = (6,7,4,2,2,1)
if type(rt) == str:
    rt = rt.split(',') 

#### Setting the Hydrogen presentation in Graphs##########   
Hydrogen_presentation = '0'
while any(i for i in range(1,4) if i == int(Hydrogen_presentation)) != True:
    Hydrogen_presentation = input('1:c_H (Default) | 2: wt% H | 3: c_H and wt% H \ninput:')    
    if not Hydrogen_presentation or Hydrogen_presentation == None:
        Hydrogen_presentation = '1'
 ############################################################   
 
 
#######################################################################################
#  Cria pasta para salvar os arquivos da liga e arquivo txt com os dados principais   #  
#######################################################################################  
if not os.path.isdir(f'./{alloy}'): #Verifica se a pasta já existe e se não existe cria ela com a proxima linha
    os.mkdir(f'./{alloy}')
    
dnow = datetime.now()        
my_path = os.path.abspath(alloy)
my_file = f'./{alloy}.txt'
if not os.path.isfile(os.path.join(my_path, my_file)):
    file = open(os.path.join(my_path, my_file), 'w+')
    file.close()
file = open(os.path.join(my_path, my_file),'a')
file.write(f"Alloy: {alloy} | {dnow.strftime('%d/%m/%Y %H:%M')} | \n")
file.write(f"hM[alpha]= {hm_alpha} | hM[beta]= {hm_beta} | hM[delta]= {hm_delta} |\n")
file.write(f"HM[alpha] = {H_alpha}| HM[beta] = {H_beta} | HM[delta] = {H_delta} |\n")
file.write(" Values for Theta and r is in the format: structure(theta,r)\n")
file.write(f" alpha{float(rt[0]),float(rt[1])}, beta{float(rt[2]),float(rt[3])}, delta{float(rt[4]),float(rt[5])}\n")
file.write(f"cH step: {cH_step}\n")
file.close()
########################################################################################



##### Begin of calculations ##################    
begin = time.time()
for i in range(1,int(2/cH_step +1)):
    C_H.append(truncate(cH_step*i,4))
for T in temperature:
    alpha_S = []
    beta_S = []
    delta_S = []
    alpha_H = []
    beta_H = []
    delta_H = []
    alpha_G = []
    beta_G = []
    delta_G = []
    alpha_mu_M = []
    beta_mu_M = []
    delta_mu_M = []
    alpha_mu_H = []
    beta_mu_H = []
    delta_mu_H = []
    for c_H in C_H:
        #c_H = 0.001 * c_H
        
        t_alpha = float(rt[0])
        r_alpha = float(rt[1])
        t_beta = float(rt[2])
        r_beta = float(rt[3])
        t_delta = float(rt[4])
        r_delta = float(rt[5])
        t=T/1000
        A = 33.066178
        B = -11.363417
        C = 11.432816
        D = -2.772874
        E = -0.158558
        G = 172.707974
        S_0 = A*ln(t) + B*t + C*(t**2)/2 + D*(t**3)/3 - E/(2*t**2) + G # Hydrogen standard entropy J/mol of H2- NIST
        S_0 = S_0/1000 #kJ/mol of H2
        #S_0 = R * 7/2 * (1 + ln(T/9.2)) # Calculation of Hydrogen standard entropy - Fukay
        
        
        if (c_H < t_alpha/r_alpha): 
            dH_m_alpha = H_alpha + c_H * hm_alpha
            dS_m_alpha = - R * (c_H * ln(c_H/(t_alpha - (c_H*(r_alpha - 1))))+ (t_alpha - r_alpha * c_H)*ln((t_alpha - r_alpha * c_H)/(t_alpha - (c_H*(r_alpha-1))))) - (c_H/2)*S_0
            dG_m_alpha = dH_m_alpha - T * dS_m_alpha
            mu_H_alpha = hm_alpha - T * (- R* ln(c_H*(t_alpha - (c_H*(r_alpha-1)))**(r_alpha-1)/(t_alpha - (c_H*r_alpha))**r_alpha)- S_0/2)
            mu_M_alpha = dG_m_alpha - c_H * mu_H_alpha
            alpha_H.append([c_H,dH_m_alpha])
            alpha_S.append([c_H,dS_m_alpha]) 
            alpha_G.append([c_H,dG_m_alpha])
            alpha_mu_H.append([c_H,mu_H_alpha])
            alpha_mu_M.append([c_H,mu_M_alpha])
            valoresalpha[c_H] = [mu_H_alpha,mu_M_alpha]
                 
            
            
        if (c_H < t_beta/r_beta):
            dH_m_beta = H_beta + c_H * hm_beta
            dS_m_beta = - R * (c_H * ln(c_H/(t_beta - (c_H*(r_beta-1))))+ (t_beta - r_beta * c_H)*ln((t_beta - r_beta * c_H)/(t_beta - (c_H*(r_beta-1))))) - (c_H/2)*S_0
            dG_m_beta = dH_m_beta - T * dS_m_beta
            mu_H_beta = hm_beta - T * (- R* ln(c_H*(t_beta - (c_H*(r_beta-1)))**(r_beta-1)/(t_beta - (c_H*r_beta))**r_beta)- S_0/2)
            mu_M_beta = dG_m_beta - c_H * mu_H_beta
            beta_H.append([c_H,dH_m_beta])
            beta_S.append([c_H,dS_m_beta])
            beta_G.append([c_H,dG_m_beta]) 
            beta_mu_H.append([c_H,mu_H_beta])
            beta_mu_M.append([c_H,mu_M_beta])
            valoresbeta[c_H] = [mu_H_beta,mu_M_beta] 
          
        
        if (c_H < t_delta/r_delta):
            dH_m_delta = H_delta + c_H * hm_delta
            dS_m_delta = - R * (c_H * ln(c_H/(t_delta - (c_H*(r_delta-1))))+ (t_delta - r_delta * c_H)*ln((t_delta - r_delta * c_H)/(t_delta - (c_H*(r_delta-1))))) - (c_H/2)*S_0
            dG_m_delta = dH_m_delta - T * dS_m_delta
            mu_H_delta =  hm_delta - T * (- R* ln(c_H*(t_delta - (c_H*(r_delta-1)))**(r_delta-1)/(t_delta - (c_H*r_delta))**r_delta)- S_0/2)
            mu_M_delta = dG_m_delta - c_H * mu_H_delta  
            delta_H.append([c_H,dH_m_delta])
            delta_S.append([c_H,dS_m_delta])
            delta_G.append([c_H,dG_m_delta])  
            delta_mu_H.append([c_H,mu_H_delta])
            delta_mu_M.append([c_H,mu_M_delta]) 
            valoresdelta[c_H] = [mu_H_delta,mu_M_delta] 
       
        
        
            
    val_S_alpha[T]= np.array(alpha_S)
    val_S_beta[T]= np.array(beta_S)
    val_S_delta[T]= np.array(delta_S)
    val_H_alpha[T]= np.array(alpha_H)
    val_H_beta[T]= np.array(beta_H)
    val_H_delta[T]= np.array(delta_H)
    val_G_alpha[T]= np.array(alpha_G)
    val_G_beta[T]= np.array(beta_G)
    val_G_delta[T]= np.array(delta_G)
    val_mu_H_alpha[T]= np.array(alpha_mu_H)
    val_mu_H_beta[T]= np.array(beta_mu_H)
    val_mu_H_delta[T]= np.array(delta_mu_H)
    val_mu_M_alpha[T]= np.array(alpha_mu_M)
    val_mu_M_beta[T]= np.array(beta_mu_M)
    val_mu_M_delta[T]= np.array(delta_mu_M)
    
    
 ############## Searching phase equilibrium ############3
    saveab = []
    savead = []
    savebd = []
    converge = 0.01
    for c_H1 in C_H:
        for c_H2 in C_H:
            if (valoresalpha.__contains__(c_H1) and valoresbeta.__contains__(c_H2) and ((c_H1 < c_H2))):
                if (truncate(abs(valoresalpha[c_H1][0] - valoresbeta[c_H2][0]),2)<=converge and truncate(abs(valoresalpha[c_H1][1] - valoresbeta[c_H2][1]),2)<=converge):
                    saveab.append([c_H1,c_H2])
            if (valoresbeta.__contains__(c_H1) and valoresdelta.__contains__(c_H2)):
                if (truncate(abs(valoresbeta[c_H1][0] - valoresdelta[c_H2][0]),2)<=converge and truncate(abs(valoresbeta[c_H1][1] - valoresdelta[c_H2][1]),2)<=converge):
                    savebd.append([c_H1,c_H2])
            if (valoresalpha.__contains__(c_H1) and valoresdelta.__contains__(c_H2) and ((c_H1 < c_H2)>0.01)):
                if (truncate(abs(valoresalpha[c_H1][0] - valoresdelta[c_H2][0]),2)<=converge and truncate(abs(valoresalpha[c_H1][1] - valoresdelta[c_H2][1]),2)<=converge):
                    savead.append([c_H1,c_H2])

    equilibrium_all = {'ab': saveab, 'ad' : savead, 'bd' : savebd}
    if not savead and any(saveab):
        equilibrium_real = 'ab'
    elif not saveab and any(savead):
        equilibrium_real = 'ad'
    elif not saveab and not savead:
        equilibrium_real = 'non equilibrium'
    else:
        if any(savebd) and saveab[0][1] > savebd[0][0]:
            equilibrium_real = 'ad'
        else:
            minimo= min(saveab[0][0],savead[0][0])
            if minimo == saveab[0][0]:
                equilibrium_real = 'ab'
            if minimo == savead[0][0]:
                equilibrium_real = 'ad'

    print(equilibrium_real)
    verify_equilibrium.append(equilibrium_real)
    p_1 = []
    p_2 = []
    p_3 = []
    p_4 = []
    p_5 = []
    P_0 = 1
    lim_delta = 1.99
    
    if H_delta <0:
        if not equilibrium_all['ab'] and not equilibrium_all['ad']:
            for c_H in C_H:
                if c_H < t_alpha/r_alpha: #intervalo fase alpha
                    mu_alpha = float(val_mu_H_alpha[T][np.where(val_mu_H_alpha[T]==c_H)[0],1])
                    lnpH = 2* mu_alpha/(R*T)
                    peq = P_0*exp(lnpH)
                    p_1.append([c_H,peq])   
                    p_2 = [[c_H,peq]]
                    p_3 = [[c_H,peq]] 
                    p_4 = [[c_H,peq]]
                    p_5 = [[c_H,peq]] 
            p_total[T] = np.array(p_1) 
            p_1_wt = p_1
            p_2_wt = p_2
            p_3_wt = p_3
            p_4_wt = p_4
            p_5_wt = p_5
            p_1 = np.array(p_1)
            H_plat_ab = val_H_alpha[T][np.where(val_H_alpha[T]==t_alpha/(2*r_alpha))[0],1]
            S_plat_ab = val_S_alpha[T][np.where(val_S_alpha[T]==t_alpha/(2*r_alpha))[0],1]
            H_plat[T] = [H_plat_ab,None]
            S_plat[T] = [S_plat_ab,None]
            p_notequilibrium.append([1/T,ln(p_1[np.where(p_1[0:,0]==t_alpha/(2*r_alpha))[0],1])])
           
           
           
    if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):
    
        if equilibrium_real == 'ad': 
            ##########Reta tangentes dG#########
            line_ad_coef = [valoresalpha[equilibrium_all[equilibrium_real][0][0]][0],valoresalpha[equilibrium_all[equilibrium_real][0][0]][1]]
            line_ad = []
            for c_H in C_H:
                if equilibrium_all[equilibrium_real][0][0] - .05 <= c_H <= equilibrium_all[equilibrium_real][0][1] + 0.05: 
                    y = line_ad_coef[0] * c_H + line_ad_coef[1]
                    x = c_H
                    line_ad.append([x,y])
            line_ad = np.array(line_ad)
            ##############Cálculo PCI################3
            for c_H in C_H:
                if c_H < equilibrium_all[equilibrium_real][0][0]: # intervalo fase alpha
                    mu_alpha = float(val_mu_H_alpha[T][np.where(val_mu_H_alpha[T]==c_H)[0],1])
                    lnpH = 2* mu_alpha/(R*T)
                    peq = P_0*exp(lnpH)
                    p_1.append([c_H,peq])
                if  equilibrium_all[equilibrium_real][0][0]<c_H<equilibrium_all[equilibrium_real][0][1]: # plateau alpha - delta
                    p_2.append([c_H,peq])
                if  equilibrium_all[equilibrium_real][0][1]<c_H< C_H[-1]: # intervalo fase delta
                    #print(c_H,val_mu_H_delta[T][np.where(val_mu_H_delta[T]==c_H)[0],1])
                    mu_delta = float(val_mu_H_delta[T][np.where(val_mu_H_delta[T]==c_H)[0],1])
                    lnpH = 2* mu_delta/(R*T)
                    peq = P_0*exp(lnpH)
                    p_3.append([c_H,peq])
                    p_4 = [[c_H,peq]]
                    p_5 = [[c_H,peq]]  
            p_total[T] = np.array(p_1 + p_2 + p_3)
            p_1_wt = p_1
            p_2_wt = p_2
            p_3_wt = p_3
            p_4_wt = p_4
            p_5_wt = p_5
            H_plat[T] = [Hplat(val_H_alpha,equilibrium_all[equilibrium_real][0][0],val_H_delta,equilibrium_all[equilibrium_real][0][1],T), None]
            S_plat[T] = [Splat(val_S_alpha,equilibrium_all[equilibrium_real][0][0],val_S_delta,equilibrium_all[equilibrium_real][0][1],T), None]
            p_plateau_ad = exp(2* ((H_plat[T][0]/(R*T) - (S_plat[T][0]/R))))
            for i in range(0,len(p_2)):
                p_2[i][1] = p_plateau_ad
            for i in range(0,len(p_1)):
                if p_1[i][1] > p_plateau_ad:
                    p_1[i][1] = p_plateau_ad
            for i in range(0,len(p_3)):
                if p_3[i][1] < p_plateau_ad:
                    p_3[i][1] = p_plateau_ad
            p_1 = np.array(p_1)
            p_2 = np.array(p_2)
            p_3 = np.array(p_3)
            p_equilibrium_ad.append([1/T,ln(p_plateau_ad)])
            p_plat[T] = [p_plateau_ad,None]
            
            
        if equilibrium_real == 'ab':
            #equilibrium_ab = equilibrium_real
            lim_ab = equilibrium_all[equilibrium_real][0][1]
            
            if not equilibrium_all['bd'] or saveab[0][1] > savebd[0][0]:
                for c_H in C_H:
                    if c_H < equilibrium_all[equilibrium_real][0][0]: #intervalo fase alpha
                        mu_alpha = float(val_mu_H_alpha[T][np.where(val_mu_H_alpha[T]==c_H)[0],1])
                        lnpH = 2* mu_alpha/(R*T)
                        peq = P_0*exp(lnpH)
                        p_1.append([c_H,peq])
                    if  equilibrium_all[equilibrium_real][0][0]<c_H<equilibrium_all[equilibrium_real][0][1]: #plateau alpha-beta
                        p_2.append([c_H,peq])
                    if lim_ab < c_H < t_beta/r_beta: #intervalo fase beta
                        mu_beta = float(val_mu_H_beta[T][np.where(val_mu_H_beta[T]==c_H)[0],1])
                        lnpH = 2* mu_beta/(R*T)
                        peq = P_0*exp(lnpH)
                        p_3.append([c_H,peq])   
                        lastcH_beta= c_H
                        if H_delta <0:
                            p_4.append([c_H,peq])
                            p_5.append([c_H,peq])
                    if t_beta/r_beta < c_H < C_H[-1]: # plateau beta-delta
                        p_4.append([c_H,peq])
                        p_5 =  [[c_H,peq]] 
                lastcH_delta = C_H[-2]
                if H_delta > 0:
                    if t_beta/r_beta == C_H[-1]: 
                            print("You must change the values of Theta and r for Beta phase")
                    H_plat_bd = Hplat(val_H_beta,lastcH_beta,val_H_delta,lastcH_delta,T)
                    S_plat_bd = Splat(val_S_beta,lastcH_beta,val_S_delta,lastcH_delta,T)
                    p_plateau_bd = exp(2* ((H_plat_bd/(R*T) - (S_plat_bd/R))))
                    p_5[0][1] = p_plateau_bd
                    for i in range(0,len(p_4)):
                        p_4[i][1] = p_plateau_bd
                    for i in range(0,len(p_3)):
                        if p_3[i][1] > p_plateau_bd:
                            p_3[i][1] = p_plateau_bd
                            
            else:
                    
                equilibrium_bd = min(equilibrium_all['bd'])
                        
                #############Retas tangentes############
                line_ab_coef = [valoresalpha[equilibrium_all[equilibrium_real][0][0]][0],valoresalpha[equilibrium_all[equilibrium_real][0][0]][1]]
                line_ab = []
                for c_H in C_H:
                    if equilibrium_all[equilibrium_real][0][0] - .05 <= c_H <= equilibrium_all[equilibrium_real][0][1] + 0.05: 
                        y = line_ab_coef[0] * c_H + line_ab_coef[1]
                        x = c_H
                        line_ab.append([x,y])
                        
                line_bd_coef = [valoresbeta[equilibrium_bd[0]][0],valoresbeta[equilibrium_bd[0]][1]]
                line_bd = []
                for c_H in C_H:
                    if equilibrium_bd[0] - .05 <= c_H <= equilibrium_bd[1] + 0.05: 
                        y = line_bd_coef[0] * c_H + line_bd_coef[1]
                        x = c_H
                        line_bd.append([x,y])       
                line_ab = np.array(line_ab)     
                line_bd = np.array(line_bd)        
                ##############Cálculos PCT#############
        
                for c_H in C_H:
                    if c_H < equilibrium_all[equilibrium_real][0][0]: #intervalo fase alpha
                        mu_alpha = float(val_mu_H_alpha[T][np.where(val_mu_H_alpha[T]==c_H)[0],1])
                        lnpH = 2* mu_alpha/(R*T)
                        peq = P_0*exp(lnpH)
                        p_1.append([c_H,peq])
                    if  equilibrium_all[equilibrium_real][0][0]<c_H<equilibrium_all[equilibrium_real][0][1]: #plateau alpha-beta
                        p_2.append([c_H,peq])
                    if lim_ab < c_H < equilibrium_bd[0]: #intervalo fase beta
                        mu_beta = float(val_mu_H_beta[T][np.where(val_mu_H_beta[T]==c_H)[0],1])
                        lnpH = 2* mu_beta/(R*T)
                        peq = P_0*exp(lnpH)
                        p_3.append([c_H,peq])
                    if equilibrium_bd[0] < c_H < equilibrium_bd[1]: # plateau beta-delta
                        p_4.append([c_H,peq])
                    if equilibrium_bd[1] <c_H< C_H[-1]: # intervalo fase delta
                        mu_delta = float(val_mu_H_delta[T][np.where(val_mu_H_delta[T]==c_H)[0],1])
                        lnpH = 2* mu_delta/(R*T)
                        peq = P_0*exp(lnpH)
                        p_5.append([c_H,peq])
                H_plat_bd = Hplat(val_H_beta,equilibrium_bd[0],val_H_delta,equilibrium_bd[1],T)
                S_plat_bd = Splat(val_S_beta,equilibrium_bd[0],val_S_delta,equilibrium_bd[1],T)
                p_plateau_bd = exp(2* ((H_plat_bd/(R*T) - (S_plat_bd/R))))
                for i in range(0,len(p_4)):
                    p_4[i][1] = p_plateau_bd
                for i in range(0,len(p_3)):
                    if p_3[i][1] > p_plateau_bd:
                        p_3[i][1] = p_plateau_bd
                for i in range(0,len(p_5)):
                    if p_5[i][1] < p_plateau_bd:
                        p_5[i][1] = p_plateau_bd
                        
            if H_delta>0:
                p_total[T] = np.array(p_1 + p_2 + p_3 + p_4 + p_5)
            elif H_delta<0:
                p_total[T] = np.array(p_1 + p_2 + p_3)
            p_1_wt = p_1
            p_2_wt = p_2
            p_3_wt = p_3
            p_4_wt = p_4
            p_5_wt = p_5
            
            H_plat_ab = Hplat(val_H_alpha,equilibrium_all[equilibrium_real][0][0],val_H_beta,equilibrium_all[equilibrium_real][0][1],T)
            S_plat_ab = Splat(val_S_alpha,equilibrium_all[equilibrium_real][0][0],val_S_beta,equilibrium_all[equilibrium_real][0][1],T)
            if not H_plat_ab == False and not S_plat_ab == False:
                p_plateau_ab = exp(2* ((H_plat_ab/(R*T) - (S_plat_ab/R))))
                for i in range(0,len(p_2)):
                    p_2[i][1] = p_plateau_ab
                for i in range(0,len(p_1)):
                    if p_1[i][1] > p_plateau_ab:
                        p_1[i][1] = p_plateau_ab
                for i in range(0,len(p_3)):
                    if p_3[i][1] < p_plateau_ab:
                        p_3[i][1] = p_plateau_ab       
            p_1 = np.array(p_1)
            p_2 = np.array(p_2)
            p_3 = np.array(p_3)
            if H_delta >0:
                p_4 = np.array(p_4)
                p_5 = np.array(p_5)
                p_equilibrium_ab.append([1/T,ln(p_plateau_ab)])
                p_equilibrium_bd.append([1/T,ln(p_plateau_bd)])
                p_plat[T] = [p_plateau_ab,p_plateau_bd]
                H_plat[T] = [H_plat_ab,H_plat_bd]
                S_plat[T] = [S_plat_ab,S_plat_bd]
            elif H_delta <0:
                if any(equilibrium_all['ab']):
                    p_plat[T] = [p_plateau_ab,None]
                    H_plat[T] = [H_plat_ab,None]
                    S_plat[T] = [S_plat_ab,None]
                else:
                    None
            
        
########################## Convert c_H to wt% H ###############################     
    
    p_1_wt = np.array(p_1_wt)
    p_2_wt = np.array(p_2_wt)
    p_3_wt = np.array(p_3_wt)
    p_4_wt = np.array(p_4_wt)
    p_5_wt = np.array(p_5_wt)
   
    if type(p_4_wt)!= type(p_1_wt) or type(p_5_wt)!= type(p_1_wt):
        p_4_wt= np.array([[0,0]])
        p_5_wt = np.array([[0,0]])    
          
    alloy_atomic_mass=sum(composition_vetor * atomic_mass[0:,1])
    
            
    for i,j,k,l,m in zip_longest(p_1_wt[0:,0],p_2_wt[0:,0],p_3_wt[0:,0],p_4_wt[0:,0],p_5_wt[0:,0]):
        if i!= None:  Hconvert(p_1_wt,alloy_atomic_mass,i)
        if j!= None:  Hconvert(p_2_wt,alloy_atomic_mass,j)
        if k!= None:  Hconvert(p_3_wt,alloy_atomic_mass,k)
        if l!= None and any(p_4_wt[0]) != [0,0]:  Hconvert(p_4_wt,alloy_atomic_mass,l)
        if m!= None and any(p_5_wt[0]) != [0,0]:  Hconvert(p_5_wt,alloy_atomic_mass,m) 
        
    
    if not equilibrium_all['ab'] and not equilibrium_all['ad']:
        p_total_wt[T] = p_1_wt
    
    if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):
        if equilibrium_real == "ab":
            if H_delta>0:
                p_total_wt[T] = np.concatenate((p_1_wt,p_2_wt,p_3_wt,p_4_wt,p_5_wt))
            elif H_delta<0:
                p_total_wt[T] = np.concatenate((p_1_wt,p_2_wt,p_3_wt))
        elif equilibrium_real == "ad":
            p_total_wt[T] = np.concatenate((p_1_wt,p_2_wt,p_3_wt))
    end = time.time()
    tempo = end - begin
    print(tempo)
    
########################## Important Infos Text ########################
    my_path = os.path.abspath(alloy)
    my_file = f'./{alloy}.txt'
    file = open(os.path.join(my_path, my_file),'a')
    file.write("------------------------------------------------------- \n")
    file.write(f"Temperature = {T - 273.15} \u00B0C | {T} K \n")   
    if equilibrium_real == 'ab':
        cHalpha = equilibrium_all[equilibrium_real][0][0]
        cHbeta1 = equilibrium_all[equilibrium_real][0][1]
        if not equilibrium_all['bd'] or saveab[0][1] > savebd[0][0]:
            cHbeta2 = lastcH_beta
            cHdelta = lastcH_delta
        else:
            cHbeta2 = equilibrium_bd[0]
            cHdelta = equilibrium_bd[1]
        wtalpha = cHalpha * 1.01 *100/ (alloy_atomic_mass+ cHalpha *1.01)
        wtbeta1 = cHbeta1 * 1.01 *100/ (alloy_atomic_mass+ cHbeta1 *1.01)
        wtbeta2 = cHbeta2 * 1.01 *100/ (alloy_atomic_mass+ cHbeta2 *1.01)
        wtdelta = cHdelta * 1.01 *100/ (alloy_atomic_mass+ cHdelta *1.01)
        if any(equilibrium_all['ab']):
            file.write("Plateau alfa-beta \n")
            file.write(f"c_H[alpha] = {cHalpha} \n" )
            file.write(f"c_H[beta] = {cHbeta1} \n")
            file.write(f"wt% H[alpha] = {wtalpha} \n" )
            file.write(f"wt% H[beta] = {wtbeta1} \n")
            file.write(f"Plateau Pressure = {p_plat[T][0]:e} atm \n")
            file.write(f"Plateau Enthalpy = {H_plat[T][0]} kJ/mol H  \n")
            file.write(f"Plateau Entropy = {S_plat[T][0]} kJ/mol H \n")
        if H_delta>0:
            file.write("Plateau beta-delta \n")     
            file.write(f"c_H[beta] ={cHbeta2} \n")
            file.write(f"c_H[delta] = {cHdelta} \n" ) 
            file.write(f"wt% H[beta] = {wtbeta2} \n" )
            file.write(f"wt% H[delta] = {wtdelta} \n")
            file.write(f"Plateau Pressure = {p_plat[T][1]:e} atm \n")
            file.write(f"Plateau Enthalpy = {H_plat[T][1]}  kJ/mol H \n")
            file.write(f"Plateau Entropy = {S_plat[T][1]} kJ/mol H \n")      
    elif equilibrium_real == 'ad':  
        cHalpha = equilibrium_all[equilibrium_real][0][0]
        cHdelta = equilibrium_all[equilibrium_real][0][1]
        wtalpha = cHalpha * 1.01 *100/ (alloy_atomic_mass+ cHalpha *1.01)
        wtdelta = cHdelta * 1.01 *100/ (alloy_atomic_mass+ cHdelta *1.01)
        file.write("Plateau alfa-delta \n")
        file.write(f"c_H[alpha] = {cHalpha} \n" )
        file.write(f"c_H[delta] = {cHdelta} \n")
        file.write(f"wt% H[alpha] = {wtalpha} \n" )
        file.write(f"wt% H[delta] = {wtdelta} \n")
        file.write(f"Plateau Pressure = {p_plat[T][0]:e} atm \n")
        file.write(f"Plateau Enthalpy = {H_plat[T][0]} kJ/mol H  \n")
        file.write(f"Plateau Entropy = {S_plat[T][0]} kJ/mol H \n")   
    elif H_delta < 0 and not equilibrium_all['ab'] and not equilibrium_all['ad']:
        file.write("Equilibrium not found")
    file.write("------------------------------------------------------- \n")
    file.close()    
        
        
############################### Plot graphs and their data txt files #######################3        
        
        
    if not os.path.isdir(f'./{alloy}/{T - 273.15}'): #Verifica se a pasta já existe e se não existe cria ela com a proxima linha
        os.mkdir(f'./{alloy}/{T - 273.15}')
    my_path = os.path.abspath(alloy+ f"/{T - 273.15}") # Figures out the absolute path for you in case your working directory moves around.
    
    
    
    
    ########## Entropy Variation #################### 
    my_file = f'{alloy}_Sall_T{T - 273.15}C.png'
    # Create Figure (empty canvas)
    fig = plt.figure()

    # Add set of axes to figure
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

    # Plot on that set of axes
    axes.plot(val_S_alpha[T][0:,0], val_S_alpha[T][0:,1], 'b', label = r'$\Delta S^{\alpha}$')
    axes.plot(val_S_beta[T][0:,0], val_S_beta[T][0:,1], 'g', label = r'$\Delta S^{\beta}$')
    axes.plot(val_S_delta[T][0:,0], val_S_delta[T][0:,1], 'r', label = r'$\Delta S^{\delta}$')
    axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
    axes.set_ylabel(r"$\Delta S$ [kJ/mol]")
    axes.set_title(f'T = {T- 273.15} \u00B0C')
    axes.legend(loc = 0)
    fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')  
    
 
    
    ########## Enthalpy Variation ####################   
    my_file = f'{alloy}_Hall_T{T - 273.15}C.png'
    # Create Figure (empty canvas)
    fig = plt.figure()

    # Add set of axes to figure
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

    # Plot on that set of axes
    axes.plot(val_H_alpha[T][0:,0], val_H_alpha[T][0:,1], 'b', label = r'$\Delta H^{\alpha}$')
    axes.plot(val_H_beta[T][0:,0], val_H_beta[T][0:,1], 'g', label = r'$\Delta H^{\beta}$')
    axes.plot(val_H_delta[T][0:,0], val_H_delta[T][0:,1], 'r', label = r'$\Delta H^{\delta}$')
    axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
    axes.set_ylabel(r"$\Delta H$ [kJ/mol]")
    axes.set_title(f'T = {T- 273.15} \u00B0C')
    axes.legend(loc = 0)
    fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight') 
    
    
    ##########Variação da Energia livre de Gibbs####################
    if equilibrium_real =="ad":
    
        my_file = f'{alloy}_Gall_T{T - 273.15}C.png'
        # Create Figure (empty canvas)
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(val_G_alpha[T][0:,0], val_G_alpha[T][0:,1], 'b', label = r'$\Delta G^{\alpha}$')
        axes.plot(val_G_beta[T][0:,0], val_G_beta[T][0:,1], 'g', label = r'$\Delta G^{\beta}$')
        axes.plot(val_G_delta[T][0:,0], val_G_delta[T][0:,1], 'r', label = r'$\Delta G^{\delta}$')
        axes.plot(line_ad[0:,0], line_ad[0:,1], 'k--', lw =0.7)
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$\Delta G$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
        
    elif equilibrium_real =="ab":
    
        my_file = f'{alloy}_Gall_T{T - 273.15}C.png'
        # Create Figure (empty canvas)
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(val_G_alpha[T][0:,0], val_G_alpha[T][0:,1], 'b', label = r'$\Delta G^{\alpha}$')
        axes.plot(val_G_beta[T][0:,0], val_G_beta[T][0:,1], 'g', label = r'$\Delta G^{\beta}$')
        axes.plot(val_G_delta[T][0:,0], val_G_delta[T][0:,1], 'r', label = r'$\Delta G^{\delta}$')
        if not equilibrium_all['bd'] or saveab[0][1] > savebd[0][0]:
            None
        else:
            axes.plot(line_ab[0:,0], line_ab[0:,1], 'k--', lw =0.7)
            axes.plot(line_bd[0:,0], line_bd[0:,1], 'k--', lw = 0.7)
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$\Delta G$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight') 
    elif equilibrium_real == "non equilibrium":
        
        my_file = f'{alloy}_Gall_T{T - 273.15}C.png'
        # Create Figure (empty canvas)
        fig = plt.figure()

        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

        # Plot on that set of axes
        axes.plot(val_G_alpha[T][0:,0], val_G_alpha[T][0:,1], 'b', label = r'$\Delta G^{\alpha}$')
        axes.plot(val_G_beta[T][0:,0], val_G_beta[T][0:,1], 'g', label = r'$\Delta G^{\beta}$')
        axes.plot(val_G_delta[T][0:,0], val_G_delta[T][0:,1], 'r', label = r'$\Delta G^{\delta}$')
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$\Delta G$ [kJ/mol]")
        axes.set_title(f'T = {T- 273.15} \u00B0C')
        axes.legend(loc = 0)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
    
    
    ##########Potencial químico do Metal####################    
    my_file = f'{alloy}_miMall_T{T - 273.15}C.png'
    # Create Figure (empty canvas)
    fig = plt.figure()

    # Add set of axes to figure
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

    # Plot on that set of axes
    axes.plot(val_mu_M_alpha[T][0:,0], val_mu_M_alpha[T][0:,1], 'b', label = r'$\mu_{M}^{\alpha}$')
    axes.plot(val_mu_M_beta[T][0:,0], val_mu_M_beta[T][0:,1], 'g', label = r'$\mu_{M}^{\beta}$')
    axes.plot(val_mu_M_delta[T][0:,0], val_mu_M_delta[T][0:,1], 'r', label = r'$\mu_{M}^{\delta}$')
    axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
    axes.set_ylabel(r"$\mu_{M}$ [kJ/mol]")
    axes.set_title(f'T = {T- 273.15} \u00B0C')
    axes.legend(loc = 0)
    fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight') 

    
    ##################Potencial químico do Hidrogênio##############       
    my_file = f'{alloy}_miHall_T{T - 273.15}C.png'
    # Create Figure (empty canvas)
    fig = plt.figure()

    # Add set of axes to figure
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)

    # Plot on that set of axes
    axes.plot(val_mu_H_alpha[T][0:,0], val_mu_H_alpha[T][0:,1], 'b', label = r'$\mu_{H}^{\alpha}$')
    axes.plot(val_mu_H_beta[T][0:,0], val_mu_H_beta[T][0:,1], 'g', label = r'$\mu_{H}^{\beta}$')
    axes.plot(val_mu_H_delta[T][0:,0], val_mu_H_delta[T][0:,1], 'r', label = r'$\mu_{H}^{\delta}$')
    axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
    axes.set_ylabel(r"$\mu_{H}$ [kJ/mol]")
    axes.set_title(f'T = {T- 273.15} \u00B0C')
    axes.legend(loc = 0)
    fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight') 
    
    
    
    ################## PCI Diagrams ##############

    if Hydrogen_presentation == '1' or Hydrogen_presentation == '3':
        if H_delta <0:
            if not equilibrium_all['ab'] and not equilibrium_all['ad']:
                my_file = f'{alloy}_PCI_T{T - 273.15}C_cH.png'
                # Create Figure (empty canvas)
                fig = plt.figure()
        
                # Add set of axes to figure
                axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
                # Plot on that set of axes
                axes.plot(p_1[0:,0], p_1[0:,1], 'b', label = r'$\alpha$')
                plt.plot([], [], ' ', label= f"Plateau = {p_2[0,1]:e} atm")
                plt.yscale("log")
                axes.set_xlabel(r'$c_H$') # Notice the use of set_ to begin methods
                axes.set_ylabel(r"$P$ [atm]")
                axes.set_title(f'T = {T- 273.15} \u00B0C')
                axes.legend(loc = 0,prop={'size': 7.5})
                fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
                
        if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):
            
            if equilibrium_real == 'ad': 
                
                
                my_file = f'{alloy}_PCI_T{T - 273.15}C_cH.png'
                # Create Figure (empty canvas)
                fig = plt.figure()
        
                # Add set of axes to figure
                axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
                # Plot on that set of axes
                axes.plot(p_1[0:,0], p_1[0:,1], 'b', label = r'$\alpha$')
                axes.plot(p_2[0:,0], p_2[0:,1], 'g', label = r'$\alpha - \delta$')
                axes.plot(p_3[0:,0], p_3[0:,1], 'r', label = r'$\delta$')
                plt.plot([], [], ' ', label= f"Plateau = {p_2[0,1]:e} atm")
                plt.yscale("log")
                axes.set_xlabel(r'$c_H$') # Notice the use of set_ to begin methods
                axes.set_ylabel(r"$P$ [atm]")
                axes.set_title(f'T = {T- 273.15} \u00B0C')
                axes.legend(loc = 0,prop={'size': 7.5})
                fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
                
                
            if equilibrium_real == 'ab': 
            
                my_file = f'{alloy}_PCI_T{T - 273.15}C_cH.png'
                 # Create Figure (empty canvas)
                fig = plt.figure()
        
                # Add set of axes to figure
                axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
                # Plot on that set of axes
                axes.plot(p_1[0:,0], p_1[0:,1], 'b', label = r'$\alpha$')
                axes.plot(p_2[0:,0], p_2[0:,1], 'g', label = r'$\alpha - \beta$')
                axes.plot(p_3[0:,0], p_3[0:,1], 'r', label = r'$\beta$')
                if H_delta >0:
                    axes.plot(p_4[0:,0], p_4[0:,1], 'y', label = r'$\beta - \delta$')
                    axes.plot(p_5[0:,0], p_5[0:,1], 'tab:orange', label = r'$\delta$')
                plt.plot([], [], ' ', label= f"Plateau 1= {p_2[0,1]:e} atm")
                if H_delta >0:
                    plt.plot([], [], ' ', label= f"Plateau 2= {p_4[0,1]:e} atm")
                plt.yscale("log")
                axes.set_xlabel(r'$c_H$') # Notice the use of set_ to begin methods
                axes.set_ylabel(r"$P$ [atm]")
                axes.set_title(f'T = {T- 273.15} \u00B0C')
                axes.legend(loc = 0,prop={'size': 7.5})
                fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')

        
    if Hydrogen_presentation == '2' or Hydrogen_presentation == '3':
        if H_delta <0:
            if not equilibrium_all['ab'] and not equilibrium_all['ad']:    
                my_file = f'{alloy}_PCI_T{T - 273.15}C_mass.png'
                # Create Figure (empty canvas)
                fig = plt.figure()
        
                # Add set of axes to figure
                axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
                # Plot on that set of axes
                axes.plot(p_1_wt[0:,0], p_1_wt[0:,1], 'b', label = r'$\alpha$')
                plt.plot([], [], ' ', label= f"Plateau = {p_2[0,1]:e} atm")
                plt.yscale("log")
                axes.set_xlabel('wt% H') # Notice the use of set_ to begin methods
                axes.set_ylabel(r"$P$ [atm]")
                axes.set_title(f'T = {T- 273.15} \u00B0C')
                axes.legend(loc = 0,prop={'size': 7.5})
                fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
        
        if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):             
            if equilibrium_real == 'ad': 
                
                my_file = f'{alloy}_PCI_T{T - 273.15}C_mass.png'
                # Create Figure (empty canvas)
                fig = plt.figure()
        
                # Add set of axes to figure
                axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
                # Plot on that set of axes
                axes.plot(p_1_wt[0:,0], p_1_wt[0:,1], 'b', label = r'$\alpha$')
                axes.plot(p_2_wt[0:,0], p_2_wt[0:,1], 'g', label = r'$\alpha - \delta$')
                axes.plot(p_3_wt[0:,0], p_3_wt[0:,1], 'r', label = r'$\delta$')
                plt.plot([], [], ' ', label= f"Plateau = {p_2[0,1]:e} atm")
                plt.yscale("log")
                axes.set_xlabel('wt% H') # Notice the use of set_ to begin methods
                axes.set_ylabel(r"$P$ [atm]")
                axes.set_title(f'T = {T- 273.15} \u00B0C')
                axes.legend(loc = 0,prop={'size': 7.5})
                fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
            
            
            if equilibrium_real == 'ab': 
                
                
                my_file = f'{alloy}_PCI_T{T - 273.15}C_mass.png'
                 # Create Figure (empty canvas)
                fig = plt.figure()
        
                # Add set of axes to figure
                axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
                # Plot on that set of axes
                axes.plot(p_1_wt[0:,0], p_1_wt[0:,1], 'b', label = r'$\alpha$')
                axes.plot(p_2_wt[0:,0], p_2_wt[0:,1], 'g', label = r'$\alpha - \beta$')
                axes.plot(p_3_wt[0:,0], p_3_wt[0:,1], 'r', label = r'$\beta$')
                if H_delta >0:
                    axes.plot(p_4_wt[0:,0], p_4_wt[0:,1], 'y', label = r'$\beta - \delta$')
                    axes.plot(p_5_wt[0:,0], p_5_wt[0:,1], 'tab:orange', label = r'$\delta$')
                plt.plot([], [], ' ', label= f"Plateau 1= {p_2[0,1]:e} atm")
                if H_delta >0:
                    plt.plot([], [], ' ', label= f"Plateau 2= {p_4[0,1]:e} atm")
                plt.yscale("log")
                axes.set_xlabel('wt% H') # Notice the use of set_ to begin methods
                axes.set_ylabel(r"$P$ [atm]")
                axes.set_title(f'T = {T- 273.15} \u00B0C')
                axes.legend(loc = 0,prop={'size': 7.5})
                fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
                
            
            
                  
    
    ###################Arquivos fonte .txt dos dados obtidos########################
    
    if not os.path.isdir(f'./{alloy}/{T - 273.15}/data'): #Verifica se a pasta já existe e se não existe cria ela com a proxima linha
        os.mkdir(f'./{alloy}/{T - 273.15}/data')
    my_path = os.path.abspath(alloy+ f"/{T - 273.15}/data") # Figures out the absolute path for you in case your working directory moves around.
    

    my_file = f'{alloy}_Sa_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_S_alpha[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "S_alpha"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'S_alpha'])
    
    my_file = f'{alloy}_Sb_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_S_beta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "S_beta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'S_beta'])
    
    my_file = f'{alloy}_Sd_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_S_delta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "S_delta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'S_delta'])
    
    my_file = f'{alloy}_Ha_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_H_alpha[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "H_alpha"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'H_alpha'])
    
    my_file = f'{alloy}_Hb_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_H_beta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "H_beta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'H_beta'])
    
    my_file = f'{alloy}_Hd_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_H_delta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "H_delta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'H_delta'])
    
    my_file = f'{alloy}_Ga_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_G_alpha[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "G_alpha"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'G_alpha'])
    
    my_file = f'{alloy}_Gb_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_G_beta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "G_beta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'G_beta'])
    
    my_file = f'{alloy}_Gd_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_G_delta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "G_delta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'G_delta'])
    
    my_file = f'{alloy}_miMa_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_mu_M_alpha[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "mi_M_alpha"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'mi_M_alpha'])
    
    my_file = f'{alloy}_miMb_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_mu_M_beta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "mi_M_beta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'mi_M_beta'])
    
    my_file = f'{alloy}_miMd_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_mu_M_delta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "mi_M_delta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'mi_M_delta'])
    
    my_file = f'{alloy}_miHa_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_mu_H_alpha[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "mi_H_alpha"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'mi_H_alpha'])
    
    my_file = f'{alloy}_miHb_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_mu_H_beta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "mi_H_beta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'mi_H_beta'])
    
    my_file = f'{alloy}_miHd_T{T - 273.15}C.txt'
    dataframe = pd.DataFrame(val_mu_H_delta[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "mi_H_delta"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'mi_H_delta'])
    
    
    my_file = f'{alloy}_PCI_T{T - 273.15}C_cH.txt'
    dataframe = pd.DataFrame(p_total[T])
    dataframe = dataframe.rename(columns = { 0 : "c_H", 1 : "P(atm)"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['c_H', 'P(atm)'])
    
    my_file = f'{alloy}_PCI_T{T - 273.15}C_mass.txt'
    dataframe = pd.DataFrame(p_total_wt[T])
    dataframe = dataframe.rename(columns = { 0 : "wt% H", 1 : "P(atm)"})
    dataframe.to_csv(os.path.join(my_path, my_file), sep = '\t', encoding = "utf-8", index = False, columns=['wt% H', 'P(atm)'])
###################Van't Hoff plot ##############################

if len(temperature)>1:    
    my_path = os.path.abspath(alloy)    

    if all( i == 'ad' for i in verify_equilibrium):
        
        p_equilibrium_ad = np.array(p_equilibrium_ad)
        
        #---------- Van't Hoff H and S in the legend
        my_file = f'{alloy}_VantHoff.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
    
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
    
        # Plot on that set of axes
        axes.plot(p_equilibrium_ad[0:,0], p_equilibrium_ad[0:,1], 'k', label = r'$\alpha - \delta$')
        for T in [temperature[0], temperature[-1]]:
                plt.plot([], [], ' ', label= r"$\Delta H_{Plat} ^{\alpha - \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(H_plat[T][0])))
                plt.plot([], [], ' ', label= r"$\Delta S_{Plat} ^{\alpha - \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(S_plat[T][0])))
        axes.set_xlabel('1/T [$K^{-1}$]') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$ln(P_{eq}(atm))$")
        axes.set_title(f"Van't Hoff {alloy}")
        axes.legend(loc = 0,prop={'size': 5})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
      
        #----------- Van't Hoff equilibrium only
        my_file = f'{alloy}_VantHoff_equilibrium.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
    
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
    
        # Plot on that set of axes
        axes.plot(p_equilibrium_ad[0:,0], p_equilibrium_ad[0:,1], 'k', label = r'$\alpha - \delta$')
        axes.set_xlabel('1/T [$K^{-1}$]') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$ln(P_{eq}(atm))$")
        axes.set_title(f"Van't Hoff {alloy}")
        axes.legend(loc = 0,prop={'size': 8})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
      
        
      
        
    elif all( i == 'ab' for i in verify_equilibrium) and H_delta>0:
            
        p_equilibrium_ab = np.array(p_equilibrium_ab)
        p_equilibrium_bd = np.array(p_equilibrium_bd)
       
        
        #---------- Van't Hoff H and S in the legend
        my_file = f'{alloy}_VantHoff.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
    
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
    
        # Plot on that set of axes
        axes.plot(p_equilibrium_ab[0:,0], p_equilibrium_ab[0:,1], 'k', label = r'$\alpha - \beta$')
        for T in [temperature[0], temperature[-1]]:
            plt.plot([], [], ' ', label= r"$\Delta H_{Plat} ^{\alpha - \beta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(H_plat[T][0])))
            plt.plot([], [], ' ', label= r"$\Delta S_{Plat} ^{\alpha - \beta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(S_plat[T][0])))
        axes.plot(p_equilibrium_bd[0:,0], p_equilibrium_bd[0:,1], 'r', label = r'$\beta - \delta$')
        for T in [temperature[0], temperature[-1]]:   
            plt.plot([], [], ' ', label= r"$\Delta H_{Plat} ^{\beta- \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(H_plat[T][1])))
            plt.plot([], [], ' ', label= r"$\Delta S_{Plat} ^{\beta- \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(S_plat[T][1])))
        axes.set_xlabel('1/T [$K^{-1}$]') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$ln(P_{eq}(atm))$")
        axes.set_title(f"Van't Hoff {alloy}")
        axes.legend(loc = 0,prop={'size': 4})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
   
    
   
        #----------- Van't Hoff equilibrium only
        my_file = f'{alloy}_VantHoff_equilibrium.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
    
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
    
        # Plot on that set of axes
        axes.plot(p_equilibrium_ab[0:,0], p_equilibrium_ab[0:,1], 'k', label = r'$\alpha - \beta$')
        axes.plot(p_equilibrium_bd[0:,0], p_equilibrium_bd[0:,1], 'r', label = r'$\beta - \delta$')
        axes.set_xlabel('1/T [$K^{-1}$]') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$ln(P_{eq}(atm))$")
        axes.set_title(f"Van't Hoff {alloy}")
        axes.legend(loc = 0,prop={'size': 8})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
   
   
    
    elif H_delta>0:
        
        if not equilibrium_all['ab'] and not equilibrium_all['bd']:
            p_notequilibrium = np.array(p_notequilibrium)
        p_equilibrium_ad = np.array(p_equilibrium_ad)
        p_equilibrium_ab = np.array(p_equilibrium_ab)
        p_equilibrium_bd = np.array(p_equilibrium_bd)
        p_equilibrium = p_equilibrium_ab + p_equilibrium_bd
       
        
        #---------- Van't Hoff H and S in the legend
        my_file = f'{alloy}_VantHoff.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
        
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
        # Plot on that set of axes
        
        if not equilibrium_all['ab'] and not equilibrium_all['bd']:
            axes.plot(p_notequilibrium[0:,0], p_notequilibrium[0:,1], 'b', label = r'$\alpha - \delta$')
        
        
        if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):
            axes.plot(p_equilibrium_ad[0:,0], p_equilibrium_ad[0:,1], 'b', label = r'$\alpha - \delta$')
            for T in temperature:
                if H_plat[T][1] == None:
                    plt.plot([], [], ' ', label= r"$\Delta H_{Plat} ^{\alpha - \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(H_plat[T][0])))
                    plt.plot([], [], ' ', label= r"$\Delta S_{Plat} ^{\alpha - \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(S_plat[T][0])))
            axes.plot(p_equilibrium_ab[0:,0], p_equilibrium_ab[0:,1], 'k', label = r'$\alpha - \beta$')
            for T in temperature:
                if not H_plat[T][1] == None:
                    plt.plot([], [], ' ', label= r"$\Delta H_{Plat} ^{\alpha - \beta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(H_plat[T][0])))
                    plt.plot([], [], ' ', label= r"$\Delta S_{Plat} ^{\alpha - \beta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(S_plat[T][0])))
            axes.plot(p_equilibrium_bd[0:,0], p_equilibrium[0:,1], 'r', label = r'$\beta - \delta$')
            for T in temperature:
                if not H_plat[T][1] == None:
                    plt.plot([], [], ' ', label= r"$\Delta H_{Plat} ^{\beta - \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(H_plat[T][1])))
                    plt.plot([], [], ' ', label= r"$\Delta S_{Plat} ^{\beta- \delta}$(%.2f °C) = %.4f kJ/mol" %(T-273.15,float(S_plat[T][1])))
        axes.set_xlabel('1/T [$K^{-1}$]') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$ln(P_{eq}(atm))$")
        axes.set_title(f"Van't Hoff {alloy}")
        axes.legend(loc = 0,prop={'size': 4})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
        
        
        my_file = f'{alloy}_VantHoff_equilibrium.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
        
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
        # Plot on that set of axes
        
        if not equilibrium_all['ab'] and not equilibrium_all['bd']:
            axes.plot(p_notequilibrium[0:,0], p_notequilibrium[0:,1], 'b', label = r'$\alpha - \delta$')
        
        
        if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):
            axes.plot(p_equilibrium_ad[0:,0], p_equilibrium_ad[0:,1], 'b', label = r'$\alpha - \delta$')  
            axes.plot(p_equilibrium_ab[0:,0], p_equilibrium_ab[0:,1], 'k', label = r'$\alpha - \beta$')
            axes.plot(p_equilibrium_bd[0:,0], p_equilibrium[0:,1], 'r', label = r'$\beta - \delta$')         
            axes.set_xlabel('1/T [$K^{-1}$]') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$ln(P_{eq}(atm))$")
        axes.set_title(f"Van't Hoff {alloy}")
        axes.legend(loc = 0,prop={'size': 8})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
        
        
        
        
    ######################PCT Diagram#######################
    if Hydrogen_presentation == '1' or Hydrogen_presentation == '3': 
        
        
        #--------- Plateu pressures in the legend-----------
        my_file = f'{alloy}_PCT_cH.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
        
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
        # Plot on that set of axes
        for T in temperature:
            axes.plot(p_total[T][0:,0], p_total[T][0:,1], label = f"{T-273.15} \u00B0C")
            if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):
                if p_plat[T][1] == None:
                    plt.plot([], [], ' ', label= f"Plateau = {p_plat[T][0]:e} atm")
                else:
                    plt.plot([], [], ' ', label= f"Plateau 1 = {p_plat[T][0]:e} atm")  
                    plt.plot([], [], ' ', label= f"Plateau 2 = {p_plat[T][1]:e} atm")
        plt.yscale("log")
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$P$ [atm]")
        axes.set_title(f'PCT {alloy}')
        axes.legend(loc = 0,prop={'size': 5})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
   
    
       #-------- Only temperatures in the legend
        my_file = f'{alloy}_PCT_cH_Tonly.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
        
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
        # Plot on that set of axes
        for T in temperature:
            axes.plot(p_total[T][0:,0], p_total[T][0:,1], label = f"{T-273.15} \u00B0C")
        current_handles, current_labels = plt.gca().get_legend_handles_labels()
        plt.yscale("log")
        axes.set_xlabel('$c_H$') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$P$ [atm]")
        axes.set_title(f'PCT {alloy}')
        a,b = search(current_labels,str(int(temperature[0]-273.15)))[0]
        current_labels.insert(0,current_labels[a])
        del(current_labels[a+1])
        current_handles.insert(0,current_handles[a])
        del(current_handles[a+1])
        axes.legend(current_handles,current_labels,loc = 0,prop={'size': 7}, ncol=2)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')
   

    if Hydrogen_presentation == '2' or Hydrogen_presentation == '3': 

        #------- Plateau pressure in the legend
        my_file = f'{alloy}_PCT_mass.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
        
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
        # Plot on that set of axes
        for T in temperature:
            axes.plot(p_total_wt[T][0:,0], p_total_wt[T][0:,1], label = f"{T-273.15} \u00B0C")
            if any(equilibrium_all['ab']) or any(equilibrium_all['ad']):    
                if p_plat[T][1] == None:
                    plt.plot([], [], ' ', label= f"Plateau = {p_plat[T][0]:e} atm")
                else:
                    plt.plot([], [], ' ', label= f"Plateau 1 = {p_plat[T][0]:e} atm")  
                    plt.plot([], [], ' ', label= f"Plateau 2 = {p_plat[T][1]:e} atm")
        plt.yscale("log")
        axes.set_xlabel('wt% H') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$P$ [atm]")
        axes.set_title(f'PCT {alloy}')
        axes.legend(loc = 0,prop={'size': 5})
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')


       #-------- Only temperatures in the legend
        my_file = f'{alloy}_PCT_mass_Tonly.png'
        # Create Figure (empty canvas)
        fig = plt.figure()
        
        # Add set of axes to figure
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
        
        # Plot on that set of axes
        for T in temperature:
            axes.plot(p_total_wt[T][0:,0], p_total_wt[T][0:,1], label = f"{T-273.15} \u00B0C")
        current_handles, current_labels = plt.gca().get_legend_handles_labels()
        plt.yscale("log")
        axes.set_xlabel('wt% H') # Notice the use of set_ to begin methods
        axes.set_ylabel(r"$P$ [atm]")
        axes.set_title(f'PCT {alloy}')
        a,b = search(current_labels, str(int(temperature[0]-273.15)))[0]
        current_labels.insert(0,current_labels[a])
        del(current_labels[a+1])
        current_handles.insert(0,current_handles[a])
        del(current_handles[a+1])
        axes.legend(current_handles,current_labels,loc = 0,prop={'size': 7}, ncol=2)
        fig.savefig(os.path.join(my_path, my_file), dpi=200, bbox_inches='tight')


end = time.time()
tempo = end - begin
print(tempo)
print("Calculate was finished")

# In[10

# In[ ]:

    