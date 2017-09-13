import numpy as np
import matplotlib.pyplot as plt
from parameters import *
import functions
from math import degrees

######################
#### Objeto TSR ######
######################

class TsrObj:
    def  __init__(self, tsr, num_naca):
        self.tsr = tsr
        self.teta = np.arange(0,180+180/float(n),180/float(n))
        self.u_linha = np.zeros(n)
        self.alfa = np.zeros(n)
        self.re = np.zeros(n)
        self.f_l = np.zeros(n)
        self.f_d = np.zeros(n)
        self.f_n = np.zeros(n)
        self.f_teta = np.zeros(n)
        self.t = np.zeros(n)
        self.f_r_mov = np.zeros(n)
        self.f_r_perp = np.zeros(n)
        self.naca = dictNACA[num_naca]

    def calcular_tsr(self):
        ############################
        #### Inicio do Calculo #####
        ############################

        #Inicializa os arrays
        ar_teta = self.teta
        ar_u_linha = np.zeros(ar_teta.size)
        ar_alfa = np.zeros(ar_teta.size)
        ar_re = np.zeros(ar_teta.size)
        ar_f_l = np.zeros(ar_teta.size)
        ar_f_d = np.zeros(ar_teta.size)
        ar_f_n = np.zeros(ar_teta.size)
        ar_f_teta = np.zeros(ar_teta.size)
        ar_t = np.zeros(ar_teta.size)
        ar_f_r_mov = np.zeros(ar_teta.size)
        ar_f_r_perp = np.zeros(ar_teta.size)

        for i in range(ar_teta.size):
            resultado = self.definir_u_linha(ar_teta[i])
            u_linha = resultado[0]
            alfa = resultado[1]
            re = resultado[2]

            ar_u_linha[i] = u_linha
            ar_alfa[i] = degrees(alfa)
            ar_re[i] = re
            if alfa >= 0 :
                ar_f_l[i] = self.interpolar_cl(alfa)*0.5*ro*a*(u_linha**2)
            else:
                ar_f_l[i] = -self.interpolar_cl(alfa)*0.5*ro*a*(u_linha**2)
            ar_f_d[i] = self.interpolar_cd(alfa)*0.5*ro*a*(u_linha**2)
            ar_f_n[i] = ar_f_d[i]*sin(alfa) + ar_f_l[i]*cos(alfa)
            ar_f_teta[i] = ar_f_d[i]*cos(alfa) + ar_f_l[i]*sin(alfa)
            ar_t[i] = ar_f_teta[i]*r/n
            ar_f_r_mov[i] = (ar_f_teta[i]*sin(ar_teta[i])+ar_f_n[i]*cos(ar_teta[i]))/n
            ar_f_r_perp[i] = (ar_f_teta[i]*cos(ar_teta[i])+ar_f_n[i]*sin(ar_teta[i]))/n
            
            print i, u_linha, alfa, re

            if i == 89:
                print i

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
        self.u_linha = ar_u_linha
        self.alfa = ar_alfa
        self.re = ar_re
        self.f_l = ar_f_l 
        self.f_d = ar_f_d
        self.f_n = ar_f_n
        self.f_teta = ar_f_teta 
        self.t = ar_t 
        self.f_r_mov = ar_f_r_mov 
        self.f_r_perp = ar_f_r_perp 
    
    def definir_u_linha (self,teta):
        u_linha = u_inf
        alfa = 0
        teta = pi*teta/180
        
        erro = 10
        count = 0
        u_linha_n_conv = 0
        while (erro > 0.001 or count < 100):
            v_teta = omega*r + u_linha*sin(teta)
            v_r = u_linha*cos(teta)
            v_res = sqrt((v_teta)**2+(v_r)**2)
            alfa = atan(v_r/v_teta)
            re = ro*v_res*c/mu
            cd = self.interpolar_cd(alfa)
            cl = self.interpolar_cl(alfa)
            if alfa >= 0:
                cx_linha = (cd*cos(alfa)-cl*sin(alfa))*sin(teta)+(cd*sin(alfa)+cl*cos(alfa))*cos(teta)
            else:
                cx_linha = (cd*cos(alfa)-cl*sin(alfa))*sin(teta)+(cd*sin(alfa)+cl*cos(alfa))*cos(teta)*-1

            #fx_linha = cx_linha*ro*(u_linha_2**2)*a/2
            u_linha_2 = u_inf/(1 + cx_linha/4)
            erro = abs((u_linha-u_linha_2)/u_linha_2)
            
            if erro > 0.001:
                u_linha = u_linha_2
            else:
                break
            
            count += 1
            
            if count == 99:
                print "Calculo de u_linha nao convergiu para (alfa, re) ", alfa, re
                return ((u_linha+u_linha_n_conv)/2, alfa, re)
        return (u_linha, alfa, re)

    def interpolar_cd(self,alfa):
        #Utilizando o valor absoluto de alfa pois o perfil e simetrico
        alfa_deg = degrees(abs(alfa))

        if (alfa_deg >= self.naca['alfa'][0] and alfa_deg < self.naca['alfa'][len(self.naca['alfa'])-1]) :
            #find index of closest smaller then...
            for i in range(1,len(self.naca['alfa'])):
                diferenca = alfa_deg - self.naca['alfa'][i]
                if diferenca > 0:
                    continue
                else:
                    #Interpolacao de cd com o valor anterior
                    return self.naca['cd'][i-1]+(self.naca['cd'][i]-self.naca['cd'][i-1])*(alfa_deg-self.naca['alfa'][i-1])/(self.naca['alfa'][i]-self.naca['alfa'][i-1])


        else:
            if (alfa_deg == self.naca['alfa'][len(self.naca['alfa'])-1]) :
                return self.naca['cd'][len(self.naca['alfa'])-1]
        
        #Exception alfa out of range...

    def interpolar_cl(self,alfa):
        #Utilizando o valor absoluto de alfa pois o perfil e simetrico
        alfa_deg = degrees(abs(alfa))

        if (alfa_deg >= self.naca['alfa'][0] and alfa_deg < self.naca['alfa'][len(self.naca['alfa'])-1]) :
            #find index of closest smaller then...
            for i in range(1,len(self.naca['alfa'])):
                diferenca = alfa_deg - self.naca['alfa'][i]
                if diferenca > 0:
                    continue
                else:
                    #Interpolacao de cl com o valor anterior
                    return self.naca['cl'][i-1]+(self.naca['cl'][i]-self.naca['cl'][i-1])*(alfa_deg-self.naca['alfa'][i-1])/(self.naca['alfa'][i]-self.naca['alfa'][i-1])


        else:
            if (alfa_deg == self.naca['alfa'][len(self.naca['alfa'])-1]) :
                return self.naca['cl'][len(self.naca['alfa'])-1]
        
        #Exception alfa out of range...
    def plotar_por_teta(self,y):
        plt.plot(self.teta, y)
        plt.grid(True)
        plt.show()