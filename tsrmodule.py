import numpy as np
import matplotlib.pyplot as plt
from parameters import *
import functions
from math import degrees
import logging

######################
#### Objeto TSR ######
######################

class TsrObj:
    def  __init__(self, tsr, num_naca):
        self.tsr = tsr
        self.teta = np.arange(0,360+360/float(n),360/float(n))
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
        self.naca = num_naca
        self.listaRe = listaRe

    def calcular_tsr(self):

        
        logging.info('Inicializacao das arrays do modulo tsr')
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

        logging.info('Inicializacao da iteracao para teta.')
        for i in range(ar_teta.size):
            if i == 181:
                print '180'
            logging.info('Inicializacao da iteracao para teta = ' + str(i))
            resultado = self.definir_u_linha(ar_teta[i])
            u_linha = resultado[0]
            alfa = resultado[1]
            re = resultado[2]

            ar_u_linha[i] = u_linha
            ar_alfa[i] = degrees(alfa)
            ar_re[i] = re
            if alfa >= 0 :
                ar_f_l[i]=self.interpolar('cl',alfa,re)*0.5*ro*a*(u_linha**2)
            else:
                ar_f_l[i] = -self.interpolar('cl',alfa,re)*0.5*ro*a*(u_linha**2)
            ar_f_d[i] = self.interpolar('cd',alfa,re)*0.5*ro*a*(u_linha**2)
            ar_f_n[i] = ar_f_d[i]*sin(alfa) + ar_f_l[i]*cos(alfa)
            ar_f_teta[i] = ar_f_d[i]*cos(alfa) + ar_f_l[i]*sin(alfa)
            ar_t[i] = ar_f_teta[i]*r/n
            ar_f_r_mov[i] = (ar_f_teta[i]*sin(ar_teta[i])+ar_f_n[i]*cos(ar_teta[i]))/n
            ar_f_r_perp[i] = (ar_f_teta[i]*cos(ar_teta[i])+ar_f_n[i]*sin(ar_teta[i]))/n
            
            logging.info('Resultado para teta = ' + str(i))
            logging.info('f_l = ' + str(ar_f_l[i]))
            logging.info('f_d = ' + str(ar_f_d[i]))
            logging.info('f_n = ' + str(ar_f_n[i]))
            logging.info('f_teta = ' + str(ar_f_teta[i]))
            logging.info('f_t = ' + str(ar_t[i]))
            logging.info('f_r_mov = ' + str(ar_f_r_mov[i]))
            logging.info('f_r_perp = ' + str(ar_f_r_perp[i]))

        logging.info('Fim da iteracao para teta.')
        t_total = np.sum(ar_t)
        p = t_total*omega
        c_p = p/p_max
        f_r_mov_total = np.sum(ar_f_r_mov)
        f_r_perp_total = np.sum(ar_f_r_perp)

        logging.info('Resultados:')
        logging.info('t_total = ' + str(t_total))
        logging.info('p = ' + str(p))
        logging.info('c_p = ' + str(c_p))
        logging.info('f_r_mov_total = ' + str(f_r_mov_total))
        logging.info('f_r_perp_total = ' + str(f_r_perp_total))

        f_r = sqrt(f_r_mov_total**2+f_r_perp_total**2)  
        c_d_rotor = f_r_mov_total/(ro*(u_inf**2)*a*0.5)
        re_rotor = ro*(u_inf**2)*(2*r)/mu

        logging.info('f_r = ' + str(f_r))
        logging.info('c_d_rotor = ' + str(c_d_rotor))
        logging.info('re_rotor = ' + str(re_rotor))

        u_2 = 2*u_linha - u_inf

        logging.info('u_2 = ' + str(u_2))

        p_3 = p_atm + ro*0.5*(u_inf**2-u_linha**2)
        p_4 = p_atm + ro*0.5*(u_2**2-u_linha**2)

        logging.info('p_3 = ' + str(p_3))
        logging.info('p_4 = ' + str(p_4))

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

        logging.info('beta = ' + str(beta))
        logging.info('gama = ' + str(gama))

        f_r_x = f_r_mov_total*cos(gama+beta/2)+f_r_perp_total*sin(gama+beta/2)
        f_r_y = f_r_mov_total*sin(gama+beta/2)+f_r_perp_total*cos(gama+beta/2)

        logging.info('f_r_x = ' + str(f_r_x))
        logging.info('f_r_y = ' + str(f_r_y))

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
        logging.info('Inicializacao da iteracao para u_linha para teta = ' + str(i))
        # Estimativa inicial para u_linha = u_inf
        u_linha = u_inf
        alfa = 0
        # Definicao de teta em radianos
        teta = pi*teta/180
        
        erro = 10
        count = 0
        u_linha_n_conv = 0
        while (erro > 0.00001 or count < 100):
            v_teta = omega*r + u_linha*sin(teta)
            v_r = u_linha*cos(teta)
            v_res = sqrt((v_teta)**2+(v_r)**2)
            alfa = atan(v_r/v_teta)
            re = ro*v_res*c/mu
            cd = self.interpolar('cd', alfa, re)
            cl = self.interpolar('cl', alfa, re)
            if alfa >= 0:
                cx_linha = (cd*cos(alfa)-cl*sin(alfa))*sin(teta)+(cd*sin(alfa)+cl*cos(alfa))*cos(teta)
            else:
                cx_linha = (cd*cos(alfa)-cl*sin(alfa))*sin(teta)+(cd*sin(alfa)+cl*cos(alfa))*cos(teta)*-1

            #fx_linha = cx_linha*ro*(u_linha_2**2)*a/2
            u_linha_2 = u_inf/(1 + cx_linha/4)
            erro = abs((u_linha-u_linha_2)/u_linha_2)

            logging.info('Calculo de U_LINHA: Count = ' + str(count) + ' , erro = ' + str(erro) + ' , u_linha = ' + str(u_linha_2))
            
            if erro > 0.00001:
                u_linha = u_linha_2
            else:
                break
            
            count += 1
            
            if count == 99:
                logging.info('Calculo de U_LINHA nao convergiu para alfa = ' + str(alfa) + 'e re = ' + str(re))
                return ((u_linha+u_linha_n_conv)/2, alfa, re)

        logging.info('Calculo de U_LINHA: u_linha = ' + str(u_linha) + ', alfa = ' + str(alfa) + ', re = ' + str(re))
        return (u_linha, alfa, re)

    def interpolar(self, c, alfa, re):

        for i in range(0,len(self.listaRe)):
                diferenca = re - self.listaRe[i]
                if diferenca > 0:
                    if i == len(self.listaRe) - 1:
                        re1 = self.listaRe[i]
                        re2 = re1
                        return self.interpolar_em_alfa(alfa,re1,c)
                    continue
                else:
                    re2 = self.listaRe[i]
                    re1 = self.listaRe[i-1]

        c_alfa_re2 = self.interpolar_em_alfa(alfa, re2, c)
        c_alfa_re1 = self.interpolar_em_alfa(alfa, re1, c)

        return c_alfa_re2*(re1-re)/(re1-re2)+c_alfa_re1*(re-re2)/(re1-re2)

    def interpolar_em_alfa(self, alfa, re, c):
        #Utilizando o valor absoluto de alfa pois o perfil e simetrico
        # dictNACA[num_naca][re]={'alfa':[], 'cd': [], 'cl':[]}
        # dictNACA[num_naca][re][alfa]
        alfa_deg = degrees(abs(alfa))
        if (alfa_deg >= dictNACA[self.naca][str(re)]['alfa'][0] and alfa_deg < dictNACA[self.naca][str(re)]['alfa'][len(dictNACA[self.naca][str(re)]['alfa'])-1]) :
            #find index of closest smaller then...
            for i in range(1,len(dictNACA[self.naca][str(re)]['alfa'])):
                diferenca = alfa_deg - dictNACA[self.naca][str(re)]['alfa'][i]
                if diferenca > 0:
                    continue
                else:
                    #Interpolacao de cd com o valor anterior
                    c1 = dictNACA[self.naca][str(re)][c][i-1]
                    c2 = dictNACA[self.naca][str(re)][c][i]
                    alfa1 = dictNACA[self.naca][str(re)]['alfa'][i-1]
                    alfa2 = dictNACA[self.naca][str(re)]['alfa'][i]

                    return c1 + (c2-c1)*(alfa_deg-alfa1)/(alfa2-alfa1)
        else:
            if (alfa_deg == dictNACA[self.naca][str(re)]['alfa'][len(dictNACA[self.naca][str(re)]['alfa'])-1]) :
                return dictNACA[self.naca][str(re)][c][len(dictNACA[self.naca][str(re)]['alfa'])-1]