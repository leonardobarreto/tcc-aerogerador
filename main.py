print "Hello World"

from math import *
from parameters import *
import tsrmodule
import functions
import time

### Inicio do Programa ####
numero_naca = '0012'
t = tsrmodule.TsrObj(6,'0012_160000')
t.calcular_tsr()
#t.plotar(t.teta,t.f_l)

tituloResultado = "Resultado " + time.strftime("%Y%m%d-%H%M%S")
ws3 = wb.create_sheet(title=tituloResultado)

ws3.cell(column=1, row=1, value="teta")
ws3.cell(column=2, row=1, value="alfa")
ws3.cell(column=3, row=1, value="re")
ws3.cell(column=4, row=1, value="u_linha")
ws3.cell(column=5, row=1, value="f_l")
ws3.cell(column=6, row=1, value="f_d")
ws3.cell(column=7, row=1, value="f_n")
ws3.cell(column=8, row=1, value="f_teta")
ws3.cell(column=9, row=1, value="t")
ws3.cell(column=10, row=1, value="f_r_mov")
ws3.cell(column=11, row=1, value="f_r_perp")

for row in range(0, t.teta.size):
    ws3.cell(column=1, row=row+2, value=t.teta[row])
    ws3.cell(column=2, row=row+2, value=t.alfa[row])
    ws3.cell(column=3, row=row+2, value=t.re[row])
    ws3.cell(column=4, row=row+2, value=t.u_linha[row])
    ws3.cell(column=5, row=row+2, value=t.f_l[row])
    ws3.cell(column=6, row=row+2, value=t.f_d[row])
    ws3.cell(column=7, row=row+2, value=t.f_n[row])
    ws3.cell(column=8, row=row+2, value=t.f_teta[row])
    ws3.cell(column=9, row=row+2, value=t.t[row])
    ws3.cell(column=10, row=row+2, value=t.f_r_mov[row])
    ws3.cell(column=11, row=row+2, value=t.f_r_mov[row])

wb.save('dadosNACA.xlsx')
print "Bye world"