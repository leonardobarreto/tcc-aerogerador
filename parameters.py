from math import *
from openpyxl import load_workbook

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



####################
### Perfil NACA ####
####################

wb = load_workbook(filename = 'dadosNACA.xlsx')
ws = wb[wb.sheetnames[0]]
ws_range = str.split(ws.dimensions,':')

#Define o Numero de valores de alfa e o numero pares (cd,cl)
alfa_range = int(ws_range[1][1:])
naca_range = ((ord(ws_range[1][:1].lower()) - 96)-1)/2
#Linha 1: Valores Naca
#Linha 2: Valores Reynolds
#Linha 3: Titulos
#Coluna 1: Valores Alfa
dictNACA = {}
for naca in range(1,naca_range+1):
    num_naca = str(ws[chr(97+2*naca).upper()+str(1)].value)
    re = str(ws[chr(97+2*naca).upper()+str(2)].value)
    naca_name = num_naca+"_"+re
    dictNACA[naca_name]={'alfa':[],'cd':[],'cl':[]}
    for i in range(4,alfa_range+1):
        dictNACA[naca_name]['alfa'].append(float(ws['A'+str(i)].value))
        dictNACA[naca_name]['cl'].append(float(ws[chr(96+2*naca).upper() + str(i)].value))
        dictNACA[naca_name]['cd'].append(float(ws[chr(96+2*naca+1).upper() + str(i)].value))

print "Hi"