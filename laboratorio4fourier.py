#Punto 1 y 3
"""

#Librerias
import matplotlib.pyplot as plt 
#Importamos todo el modulo sympy
from sympy import *
#ademas importamos las variables simbolicas t'
from sympy.abc import t
import numpy as np 
import matplotlib.pyplot as plt 
from scipy import signal 
from sympy import *
from scipy.io.wavfile import read 
import random

"""##Tren de Pulsos (POLAR )"""

#calculo de coeficinetes

#P= 0.00001    para F de 100 KHz
#P= 0.000001   Para F de 1   Mhz
#P= 0.00000001 Para F de 100 Mhz
periodo = 0.00001
frecuency = 1/periodo +0.0000001
pi=np.pi   #defino pi

w=2*pi*frecuency # defino la frecuencia angular
T2=periodo/2

# calculo de  a0
a01 = (2/periodo)*integrate(1, (t, 0, T2))
a02 = (2/periodo)*integrate(-1, (t,  T2, periodo))
a0  = a01+a02
F0=a0/2
Pot0= F0**2

pott=abs(Pot0) # potencia media

# Calculo de potencia media de la señal
P= (1/periodo)*(integrate(abs((1))**2, (t, 0, T2))) + (1/periodo)* (integrate(abs((-1))**2, (t, T2, periodo)))
P9=P*0.9

pdb=10*np.log10(float(P/0.001))   #Potencia Total en dBm
print('Potencia Media = {} -- Potencia dBm: {} \n90% de Potencia = {}\n'. format(P,pdb,P9))

print ('Coeficientes : ')
print ('\na0/2=  {} -- Serie Trigonometrica \nC0= {}----Potencia C0 : {}\n'.format(F0,F0,Pot0))

# numero de coeficientes  g 
g= 1 + 5   # 5 armonicos + DC
# vectores para almacenar los coeficientes
a=np.zeros(g)
b=np.zeros(g)
c=np.zeros(g)
c[0]=F0
# calculo de coeficientes para K diferente de 0
for k in range(1, g):
  
  a[k]= (2/periodo)*(integrate((1)  * cos( w* k * t), (t, 0, T2))) + (2/periodo)*(integrate((-1)  * cos( w* k * t), (t, T2, periodo)))
  b[k]= (2/periodo)*(integrate((1)  * sin( w* k * t), (t, 0, T2))) + (2/periodo)*(integrate((-1)  * sin( w* k * t), (t, T2, periodo)))
  c[k] = (sqrt(((a[k])**2)+((b[k])**2)))/2
  CnPot= c[k]**2  
  pott=pott + abs(2*(c[k]**2)) # sumatoria de potencia
  pow=CnPot*2                  # potencia del armónico
  fre=frecuency*k
  print ("a"+str(k)+"="+ str(a[k])+"   b"+str(k)+ "="+str(b[k])+" ----  De la serie Trigonométrica")
  print('C{} = {} --- Potencia 2*|Ck|^2 : {} ---- Frecuencia: {} HZ  -- P acumulada  = {}\n'.format(k,c[k],pow,fre,pott))

#probamos reconstruir la señal, esto NO LO PIDE EL INFORME, 
#  PARA verificar si los coeficientes calcuados son correctos

t = np.arange(0, 3*periodo, periodo/100)
w=2*pi*frecuency
w=float(w)
F1=0
for k in range(1,g):
  F1=F1 + a[k]*np.cos(w*k*t)+ b[k]*np.sin(w*(k)*t) # calcular la señal por medio de coeficientes

F1 = F0 + F1
plt.title('Tren de pulsos RECONSTRUIDO')
plt.plot(t,F1)


#Graficar Espectro de Amplitud
 
plt.figure()
plt.grid()
plt.title('Espectro de Amplitud |Ck|')
freq=np.arange(-5,8,1)
cn=np.zeros(13)
for k in range(0,11,1):  
  if(k<5):
    cn[k]=c[5-k]
  if(k==5):
    cn[k]=c[0]
  if(k>5):
    cn[k]=c[k-5]
plt.stem(freq, cn)
plt.xlim(-6,7)
plt.xlabel('K')

"""## Diente de Sierra"""

#calculo de coeficientes
#Importamos todo el modulo sympy
#P= 0.00001    para F de 100 KHz
#P= 0.000001   Para F de 1   Mhz
#P= 0.00000001 Para F de 100 Mhz
T0=0.000001    # Periodo
dutty=1
periodo=T0
frecuency=1/T0 +0.000001   #Frecuencia
pi=np.pi                   #PI 
w=2 *pi*frecuency          #Frecuencia angular
#calculo de a0, C0  
a0 = (2/periodo)*integrate((1/T0)*t, (t, 0, periodo))

F0=a0/2
Pot0= F0**2
pott=abs(Pot0) # potencia media

# Calculo de potencia media de la señal
P= (1/periodo)*(integrate(abs((1/T0)*t)**2, (t, 0, T0))) 
P9=P*0.9
pdb=10*np.log10(float(P/0.001))   #Potencia Total en dBm
print('Potencia Media = {} -- Potencia dBm: {} \n90% de Potencia = {}\n'. format(P,pdb,P9))

print ('Coeficientes : \n')
print ('a0 = {}\na0/2 => C0= {}\n'.format(a0,F0))

# numero de coeficientes  g 
g= 1 + 5   # 5 armonicos + DC
#vectores para almacenar los coeficientes
a=np.zeros(g)
b=np.zeros(g)
c=np.zeros(g)
c[0]=F0
#calculo de coeficientes para k diferente de 0
for k in range(1, g):
  
  a[k]= (2/periodo)*integrate((1/T0)*t * cos( w* k * t), (t, 0, T0))
  b[k]= (2/periodo)*integrate((1/T0)*t * sin( w* k * t), (t, 0, T0))
  c[k] = (sqrt(((a[k])**2)+((b[k])**2)))/2
  CnPot= c[k]**2
  pott=pott + abs(2*(c[k]**2))  #Sumatoria de potencia
  pow=CnPot*2                   #potencia de cada armónico +-K
  fre=frecuency*k 
  print ("a"+str(k)+"="+ str(a[k])+"   b"+str(k)+ "="+str(b[k])+" ----  De la serie Trigonométrica")
  print('C{} = {} --- Potencia 2*|Ck|^2 : {} ---- Frecuencia: {} HZ  -- P acumulada  = {}\n'.format(k,c[k],pow,fre,pott))

#probamos reconstruir la señal,  
#  PARA verificar si los coeficientes calcuados son correctos

t = np.arange(0, 3*periodo, periodo/100)
w=2*pi*frecuency
w=float(w)
F1=0
for k in range(1,g):
  F1=F1 + a[k]*np.cos(w*k*t)+ b[k]*np.sin(w*(k)*t) #sumatoria de señales cosenoidales ponderadas por los coeficinets de las serie trigonométrica

F1 = F0 + F1
plt.title('Diente de Sierra RECONSTRUIDO')
plt.plot(t,F1)


 #Graficar Espectro de Amplitud
 
plt.figure()
plt.grid()
plt.title('Espectro de Amplitud')
freq=np.arange(-5,8,1)
cn=np.zeros(13)
for k in range(0,11,1):
  if(k<5):
    cn[k]=c[5-k]
  if(k==5):
    cn[k]=c[0]
  if(k>5):
    cn[k]=c[k-5]
plt.stem(freq, cn)
plt.xlim(-6,7)
plt.xlabel('K')

"""##Triangular"""

#P= 0.00001    para F de 100 KHz
#P= 0.000001   Para F de 1   Mhz
#P= 0.00000001 Para F de 100 Mhz

periodo=0.00001  #Periodo
p=periodo
t2=periodo/2
frec = 1/periodo  #frecuencia
pi=np.pi          #PI
w=2*pi*frec       #Frecuencia angular
# calculo de a0
a01 = (2/periodo)*integrate((2/p)*t, (t, 0, t2 ))
a02 = (2/periodo)*integrate(((-2/p)*(t-p)), (t,  t2, periodo))
a0  = a01+a02
Pot0= F0**2  #Potencia en K=0
#valor dc
F0=a0/2
Pot0= F0**2 #pot del DC
pott=abs(Pot0) # potencia media acumulada

# Calculo de potencia media de la señal
P= (1/periodo)*(integrate(abs((2/p)*(t))**2, (t, 0, t2))) + (1/periodo)* (integrate(abs(((-2/p)*(t-p)))**2, (t, t2, periodo)))
P9=P*0.9

pdb=10*np.log10(float(P/0.001))   #Potencia Total en dBm
print('Potencia Media = {} -- Potencia dBm: {} \n90% de Potencia = {}\n'. format(P,pdb,P9))
print ('Coeficientes : \n')
print ('\nC0=  {} -- Serie Trigonometrica \nC0= {}----Potencia C0 : {}\n'.format(F0,F0,Pot0))


# numero de coeficientes  g 
g= 1 + 5   # 5 armónicos + DC
#vectores para almacenar los coeficientes
a=np.zeros(g)
b=np.zeros(g)
c=np.zeros(g)
c[0]=F0
for k in range(1, g):
  
  a[k]= (2/periodo)*(integrate((2/p)*(t)  * cos( w* k * t), (t, 0, t2))) + (2/periodo)*(integrate(((-2/p)*(t-p))  * cos( w* k * t), (t, t2, periodo)))
  b[k]= (2/periodo)*(integrate((2/p)*(t)  * sin( w* k * t), (t, 0, t2))) + (2/periodo)*(integrate(((-2/p)*(t-p))  * sin( w* k * t), (t, t2, periodo)))
  c[k] = (sqrt(((a[k])**2)+((b[k])**2)))/2
  CnPot= c[k]**2
  pott=pott + abs(2*(c[k]**2))
  pow=CnPot*2
  fre=frecuency*k
  print ("a"+str(k)+"="+ str(a[k])+"   b"+str(k)+ "="+str(b[k])+" ----  De la serie Trigonométrica")
  print('C{} = {} --- Potencia 2*|Ck|^2 : {} ---- Frecuencia: {} HZ  -- P acumulada  = {}\n'.format(k,c[k],pow,fre,pott))

#probamos reconstruir la señal, esto NO LO PIDE EL INFORME, 
#  PARA verificar si los coeficientes calcuados son correctos

t = np.arange(0, 3*periodo, periodo/100)
w=2*pi*frec
w=float(w)
F1=0
for k in range(1,g):
  F1=F1 + a[k]*np.cos(w*k*t)+ b[k]*np.sin(w*(k)*t) 

F1 = F0 + F1
plt.title('Triangular RECONSTRUIDa')
plt.plot(t,F1)


 #Graficar Espectro de Amplitud
 
plt.figure()
plt.grid()
plt.title('Espectro de Amplitud')
freq=np.arange(-5,8,1)
cn=np.zeros(13)
for k in range(0,11,1):
  if(k<5):
    cn[k]=c[5-k]
  if(k==5):
    cn[k]=c[0]
  if(k>5):
    cn[k]=c[k-5]
plt.stem(freq, cn)
plt.xlim(-6,7)
plt.xlabel('K')


#Punto 4
"""

#audio

import numpy as np
import scipy.fftpack as fourier
import matplotlib.pyplot as plt
import scipy.io.wavfile as waves

from scipy import signal

from scipy.fftpack import fft, fftfreq

# archivo = input('archivo de sonido:' )
archivo = 'audio.wav'
muestreo, sonido = waves.read(archivo)

#normalizando la señal 
sonido=sonido/np.max(sonido)
sonido=sonido-np.mean(sonido)

# canales: monofónico o estéreo
tamano = np.shape(sonido)
muestras = tamano[0]
m = len(tamano)
canales = 1  # monofónico
if (m>1):  # estéreo
    canales = tamano[1]
# experimento con un canal
if (canales>1):
    canal = 0
    uncanal = sonido[:,canal] 
else:
    uncanal = sonido
    

# rango de observación en segundos segmento de la señal 
inicia = 0.5
termina = 0.70
ini = 4.0
ter = 4.3
# observación en número de muestra
a = int(inicia*muestreo)
b = int(termina*muestreo)

c = int(ini*muestreo)
d = int(ter*muestreo)


#vector segmento sordo 
par= uncanal[c:d]

#vector de segmento de laseñal 
parte = uncanal[a:b]
dt = 1/muestreo
t=np.arange(0,len(sonido)/muestreo,dt)

# Salida # Archivo de audio.wav
waves.write('parte01.wav', muestreo, parte)

# Gráfica
#vector de tiempo 

t=np.arange(0,len(sonido)/muestreo,dt)
plt.plot(t,sonido)
plt.xlabel('tiempo (s)')
plt.ylabel('Amplitud')
plt.title('SEÑAL DE AUDIO',fontsize=15)
plt.show()

#transformada de fourier 
yw=fft(sonido)
freq=np.fft.fftfreq(sonido.size,d=dt)
plt.title('Transformada de fourier ',fontsize=15)
plt.xlabel('F')
plt.ylabel('Amplitud')
plt.axis(xmin=-4000,xmax=4000)
plt.plot(freq,yw)

'''
N=1024
w,Am = signal.freqz(sonido,whole=True, worN=N) # Respuesta en frecuencia de la señal original

w_transform=(w-np.pi)*muestreo/(2*np.pi)

Am_transform=(np.abs(np.fft.fftshift(Am)))**2

Am_transform=Am_transform/np.max(Am_transform)

plt.figure(figsize=(15,7))
plt.grid()
plt.plot(w_transform,Am_transform) # Grafica respuesta en frecuencia
plt.axis(xmin=-7500,xmax=7500)
plt.ylabel('Amplitud',fontsize=18) # Etiqueta eje X
plt.title('Respuesta en frecuencia de la señal de audio',fontsize=15)
plt.show()
'''

"""#Punto 5:

"""

#5 muestra de 70ms

yw=fft(parte)
freq=np.fft.fftfreq(parte .size,d=dt)
plt.title('Transformada de fourier segmento ',fontsize=15)
plt.xlabel('Frecuencia')
plt.ylabel('Amplitud')
plt.axis(xmin=-4000,xmax=4000)
plt.plot(freq,yw)

"""# Punto 8
En este caso se implementa un tren de pulsos unipolar de amplitud 1 , al cual se le variará el Dutty
"""

# Variando Dutty en Tren de Pulsos

#P= 0.00001    para F de 100 KHz
#P= 0.000001   Para F de 1   Mhz
#P= 0.00000001 Para F de 100 Mhz
periodo = 0.000001
frecuency = 1/periodo +0.0000001
T=periodo
pi=np.pi

# definimos el DUTTY deseado
D=0.50
#defino la frecuencia angular
w=2*pi*frecuency

# calculo de  a0
a01 = (2/periodo)*integrate(1, (t, 0, T*D))
#a02 = (2/periodo)*integrate(-1, (t,  T*D, periodo))
a0  = a01 #+a02
F0=a0/2
Pot0= F0**2

pott=abs(Pot0) # potencia media

# Calculo de potencia media de la señal
P= (1/periodo)*(integrate(abs((1))**2, (t, 0, T*D))) #+ (0/periodo)* (integrate(abs((-1))**2, (t, T2, periodo)))
P9=P*0.9
#pdb=10*np.log10(P)   #Potencia Total en dB
print('Potencia Media = {}\n90% de Potencia = {} \n'. format(P,P9))

print ('Coeficientes : ')
print ('\na0/2=  {} -- Serie Trigonometrica \nC0= {}----Potencia C0 : {}\n'.format(F0,F0,Pot0))

# numero de coeficientes  g 
g= 1 + 5   # 5 armonicos + DC

a=np.zeros(g)
b=np.zeros(g)
c=np.zeros(g)
c[0]=abs(F0)

for k in range(1, g):
  
  a[k]= (2/periodo)*(integrate((1)  * cos( w* k * t), (t, 0, T*D))) #+ (0/periodo)*(integrate((-1)  * cos( w* k * t), (t, T*D, periodo)))
  b[k]= (2/periodo)*(integrate((1)  * sin( w* k * t), (t, 0, T*D))) #+ (0/periodo)*(integrate((-1)  * sin( w* k * t), (t, T*D, periodo)))
  c[k] = (sqrt(((a[k])**2)+((b[k])**2)))/2
  CnPot= c[k]**2
  pott=pott + abs(2*(c[k]**2))
  pow=CnPot*2
  fre=frecuency*k
  print ("a"+str(k)+"="+ str(a[k])+"   b"+str(k)+ "="+str(b[k])+" ----  De la serie Trigonométrica")
  print('C{} = {} --- Potencia 2*|Ck|^2 : {} ---- Frecuencia: {} HZ  -- P acumulada  = {}\n'.format(k,c[k],pow,fre,pott))

#probamos reconstruir la señal 
#  PARA verificar si los coeficientes calcuados son correctos

t = np.arange(0, 3*periodo, periodo/100)
w=2*pi*frecuency
w=float(w)
F1=0
for k in range(1,g):
  F1=F1 + a[k]*np.cos(w*k*t)+ b[k]*np.sin(w*(k)*t) 

F1 = F0 + F1
plt.title("Tren de pulsos RECONSTRUIDO Dutty = "+str(D*100)+"%")
plt.plot(t,F1)


#Graficar Espectro de Amplitud
 
plt.figure()
plt.grid()
plt.title("Espectro de Amplitud |Ck| Dutty = "+str(D*100)+"%")
freq=np.arange(-5,8,1)
cn=np.zeros(13)
for k in range(0,11,1):
  if(k<5):
    cn[k]=c[5-k]
  if(k==5):
    cn[k]=c[0]
  if(k>5):
    cn[k]=c[k-5]
plt.stem(freq, cn)
plt.xlim(-6,7)
plt.xlabel('K')


