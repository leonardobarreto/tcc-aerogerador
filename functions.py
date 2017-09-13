from parameters import *
from math import *

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
        elif count == 100:
            print "O calculo de u_linha nao convergiu..."
            break
        else:
            break
        
        count += 1
    return (u_linha, alfa)