#######################################
#######################################
##  Simulacao de um Aerogerador Vertical
## Leonardo Barreto - 08/2017
#######ää##############################
#######################################

from math import *
import numpy as np
import matplotlib.pyplot as plt

#############################
##### Parametros Globais ####
#############################

#Velocidade do vento [m/s]
u_inf = 10
#TSR
tsr = 6
#Raio do aerogerador [m]
r = 0.64
#Velocidade de rotacao 
omega = tsr*u_inf/r 
#Densidade do ar [kg/m3]
ro = 1.19
#Viscodidade do ar [kg/ms] 
mu = 0.0000184
#Area
a = 1
#Corda [m]
c = 0.032
#Cd0
cd_0 = 1
#Cl0
cl_0 = 1
#Epsilon
epsilon = 0.01
#Potencia Maxima
p_max = 0.5*ro*a*(u_inf**3)
#Pressao Atmosferica [Pa]
p_atm = 101325

#####################
#### Inputs     #####
#####################
n = 180
dteta = 2*pi/n

#####################
#### Restricoes #####
#####################

re_min = 10000
re_max = 10000000

######################
#### Objeto TSR ######
######################

class TsrObj:
    def  __init__(self, tsr):
        self.tsr = tsr
        self.teta = np.arange(0,180+180/n,180/n)
        self.f_l = np.zeros(n)
        self.f_d = np.zeros(n)
        self.f_n = np.zeros(n)
        self.f_teta = np.zeros(n)
        self.t = np.zeros(n)
        self.f_r_mov = np.zeros(n)
        self.f_r_perp = np.zeros(n)

    def calcular_tsr(self):
        ############################
        #### Inicio do Calculo #####
        ############################

        #Inicializa os arrays
        ar_teta = np.arange(0,90,10)
        ar_f_l = np.zeros(ar_teta.size)
        ar_f_d = np.zeros(ar_teta.size)
        ar_f_n = np.zeros(ar_teta.size)
        ar_f_teta = np.zeros(ar_teta.size)
        ar_t = np.zeros(ar_teta.size)
        ar_f_r_mov = np.zeros(ar_teta.size)
        ar_f_r_perp = np.zeros(ar_teta.size)

        for i in range(ar_teta.size):
            resultado = definir_u_linha(ar_teta[i])
            u_linha = resultado[0]
            alfa = resultado[1]

            ar_f_l[i] = cl_0*0.5*ro*a*(u_linha**2)
            ar_f_d[i] = cd_0*0.5*ro*a*(u_linha**2)
            ar_f_n[i] = ar_f_d[i]*sin(alfa) + ar_f_l[i]*cos(alfa)
            ar_f_teta[i] = ar_f_d[i]*cos(alfa) + ar_f_l[i]*sin(alfa)
            ar_t[i] = ar_f_teta[i]*r/n
            ar_f_r_mov[i] = (ar_f_teta[i]*sin(ar_teta[i])+ar_f_n[i]*cos(ar_teta[i]))/n
            ar_f_r_perp[i] = (ar_f_teta[i]*cos(ar_teta[i])+ar_f_n[i]*sin(ar_teta[i]))/n

        t_total = np.sum(ar_t)
        p = t_total*omega
        c_p = p/p_max
        f_r_mov_total = np.sum(ar_f_r_mov)
        f_r_perp_total = np.sum(ar_f_r_perp)

        f_r = sqrt(f_r_mov_total**2+f_r_perp_total**2)  
        c_d_rotor = f_r_mov_total/(ro*(u_inf**2)*a*0.5)
        re_rotor = ro*(u_inf**2)*(2*r)/mu

        u_2 = 2*u_linha - u_inf

        p_3 = p_atm + ro*0.5*(u_inf**2-u_linha**2)
        p_4 = p_atm + ro*0.5*(u_2**2-u_linha**2)

        erro = 10
        count = 0
        beta_linha = 10
        while (erro > 0.01 and count < 100):
            k1 = (f_r**2+(p_4*a)**2+(p_3*a)**2)/(2*p_3*a*(f_r-(p_4*a*cos(beta_linha)/cos(beta_linha/2)))+2*p_4*a*f_r)
            #excecao para se k1 estiver fora do dominio
            beta = 2*np.arccos(k1)

            erro = abs((beta-beta_linha)/beta_linha)
            if erro > 0.1:
                beta_linha = beta
            else:
                break

            count +=1

        k2 = (-f_r*cos(beta/2)+p_4*a*cos(beta/2)-p_3*a)/(p_4*a*sin(beta)-f_r*(beta/2))
        #excessao para se k2 estiver fora do dominio
        gama = np.arctan(k2)

        f_r_x = f_r_mov_total*cos(gama+beta/2)+f_r_perp_total*sin(gama+beta/2)
        f_r_y = f_r_mov_total*sin(gama+beta/2)+f_r_perp_total*cos(gama+beta/2)

        ######## Fim do Calculo

        self.teta = ar_teta 
        self.f_l = ar_f_l 
        self.f_d = ar_f_d
        self.f_n = ar_f_n
        self.f_teta = ar_f_teta 
        self.t = ar_t 
        self.f_r_mov = ar_f_r_mov 
        self.f_r_perp = ar_f_r_perp 


#####################
#### Funcoes    #####
#####################
#Funcao iterativa que define a velocidade u_linha na regiao de entrada do volume de controle
#Entrada com teta em graus
def definir_u_linha (teta):
    u_linha = u_inf
    alfa = 0
    teta = pi*teta/180
    
    erro = 10
    count = 0
    while (erro > 0.001 or count < 100):
        v_teta = omega*r + u_linha*sin(teta)
        v_r = u_linha*cos(teta)
        v_res = sqrt((v_teta)**2+(v_r)**2)
        alfa = atan(v_r/v_teta)
        re = ro*v_res*c/mu
        cx_linha = (cd_0*cos(alfa)-cl_0*sin(alfa))*sin(teta)+(cd_0*sin(alfa)+cl_0*cos(alfa))*cos(teta)
        #fx_linha = cx_linha*ro*(u_linha_2**2)*a/2
        u_linha_2 = u_inf/(1 + cx_linha/4)
        erro = abs((u_linha-u_linha_2)/u_linha_2)
        
        if erro > 0.001:
            u_linha = u_linha_2
        else:
            break
        
        print count, erro
        count += 1
    return (u_linha, alfa)

### Inicio do Programa ####

t = TsrObj(6)
t.calcular_tsr()
print t

''' ############################
#### Inicio do Calculo #####
############################

#Inicializa os arrays
ar_teta = np.arange(0,90,10)
ar_f_l = np.zeros(ar_teta.size)
ar_f_d = np.zeros(ar_teta.size)
ar_f_n = np.zeros(ar_teta.size)
ar_f_teta = np.zeros(ar_teta.size)
ar_t = np.zeros(ar_teta.size)
ar_f_r_mov = np.zeros(ar_teta.size)
ar_f_r_perp = np.zeros(ar_teta.size)

for i in range(ar_teta.size):
    u_linha = u_inf
    alfa = 0
    definir_u_linha(ar_teta[i])

    ar_f_l[i] = cl_0*0.5*ro*a*(u_linha**2)
    ar_f_d[i] = cd_0*0.5*ro*a*(u_linha**2)
    ar_f_n[i] = ar_f_d[i]*sin(alfa) + ar_f_l[i]*cos(alfa)
    ar_f_teta[i] = ar_f_d[i]*cos(alfa) + ar_f_l[i]*sin(alfa)
    ar_t[i] = ar_f_teta[i]*r/n
    ar_f_r_mov[i] = (ar_f_teta[i]*sin(ar_teta[i])+ar_f_n[i]*cos(ar_teta[i]))/n
    ar_f_r_perp[i] = (ar_f_teta[i]*cos(ar_teta[i])+ar_f_n[i]*sin(ar_teta[i]))/n

t_total = np.sum(ar_t)
p = t_total*omega
c_p = p/p_max
f_r_mov_total = np.sum(ar_f_r_mov)
f_r_perp_total = np.sum(ar_f_r_perp)

f_r = sqrt(f_r_mov_total**2+f_r_perp_total**2)  
c_d_rotor = f_r_mov_total/(ro*(u_inf**2)*a*0.5)
re_rotor = ro*(u_inf**2)*(2*r)/mu

u_2 = 2*u_linha - u_inf

p_3 = p_atm + ro*0.5*(u_inf**2-u_linha**2)
p_4 = p_atm + ro*0.5*(u_2**2-u_linha**2)

erro = 10
count = 0
beta_linha = 10
while (erro > 0.01 and count < 100):
    k1 = (f_r**2+(p_4*a)**2+(p_3*a)**2)/(2*p_3*a*(f_r-(p_4*a*cos(beta_linha)/cos(beta_linha/2)))+2*p_4*a*f_r)
    #excecao para se k1 estiver fora do dominio
    beta = 2*np.arccos(k1)

    erro = abs((beta-beta_linha)/beta_linha)
    if erro > 0.1:
        beta_linha = beta
    else:
        break

    count +=1

k2 = (-f_r*cos(beta/2)+p_4*a*cos(beta/2)-p_3*a)/(p_4*a*sin(beta)-f_r*(beta/2))
#excessao para se k2 estiver fora do dominio
gama = np.arctan(k2)

f_r_x = f_r_mov_total*cos(gama+beta/2)+f_r_perp_total*sin(gama+beta/2)
f_r_y = f_r_mov_total*sin(gama+beta/2)+f_r_perp_total*cos(gama+beta/2) '''