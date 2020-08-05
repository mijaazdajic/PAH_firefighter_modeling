# -*- coding: utf-8 -*-
"""Final Project.py

Modelling skin absorption using the Euler method vs. equation from Jepps et al. (2013)
"""

import numpy as np
import matplotlib.pyplot as plt
from google.colab import files


# time period of excretion = 43 min fire + 18 hours afterwards = 1123 mins

tmax=1000           #hours
dt = 50              #time step 
t=np.arange(0,tmax,dt) 

############  Constants ############################################## 
Ksc = 0.177         #Partioning coefficient 
D_sc= 0.000093      #Diffusion coefficient in cm^2/h
Area = 1            #total skin surface area in cm^2
hsc = 0.001485      #skin depth in cm (14.85 um)
Cskin = 2.36        #max dose in ng/cm2
#########################################################################

Skin = np.zeros(len(t))
nap_jepps = np.zeros(len(t))
Skin[0]=Cskin  #exposure given in a single dose
nap_jepps[0]=Cskin

############ define function to solve model using Euler #################

def ODE_dSkin(Ksc, D_sc, Area, Skin, hsc, dt): #removing t-hsc^2/6D_sc

#dS=(-A)*dt
  dSkin = - ((((Ksc*D_sc*Area*Skin)/hsc))) * dt   

  return dSkin     

############ end of function ############################################## 

for i in range(len(t)-1):   
      
      dSkin = ODE_dSkin(Ksc, D_sc, Area, Skin[i], hsc, dt)  
      Skin[i+1]= Skin[i] + dSkin

      skin_jepps=((Ksc*D_sc*Area*Skin)/hsc)*(i-((hsc**2)/(6*D_sc)))  #equation from Jepps et al. (2013)         

# plot solution
plt.plot(t,Skin, color = "black") 
plt.plot(t,skin_jepps, 'red') 
plt.yscale('log')
#plt.xscale('log')
plt.xlabel("Time, hours")
plt.ylabel('Naphthalene, ng') 
plt.legend(['Skin surface, present study','Skin surface, Jepps et al. (2013)'])


plt.savefig("abc.svg", dpi=300, bbox_inches='tight')
files.download("abc.svg")

"""Modelling naphthalene exposure with respect to time, this is the full model"""

import numpy as np
import matplotlib.pyplot as plt
from google.colab import files


# time period of absorption = 0 - 180 mins (approx time until shower)
# time period of excretion = 43 min fire + 18 hours afterwards = 1123 mins

tmax=1123            #minutes
dt = 20               #time step 
t=np.arange(0,tmax,dt) 

############  Constants ############################################## 
Ksc = 0.177          #partioning coefficient
D_sc= 0.000093/60    #Diffusion coefficient in cm^2/h, , converted to cm^2/min 
Area = 1             #total skin surface area in cm^2
hsc = 0.001485       #skin depth in cm (14.85 um)
kelim1nap= 0.201/60  #urinary elimination coefficient for 1-hydroxynaphthalene in 1/h -> divided by 60=1/min
kelim2nap= 0.288/60  #urinary elimination coefficient for 2-hydroxynaphthalene in 1/h -> divided by 60=1/min
Cskin = 2.36         #max dose in ng/cm2
#########################################################################

S = np.zeros(len(t))
Skin = np.zeros(len(t))
SC = np.zeros(len(t))
U = np.zeros(len(t))

Skin[0]=Cskin  #exposure given in a single dose
SC[0]=0
U[0]=0

############ define function to solve model using Euler #################

def ODE_dI(t, Ksc, D_sc, Area, Skin, hsc, S, kelim1nap, kelim2nap, SC,dt):

#dS=(-A)*dt
  dSkin = - ((((Ksc*D_sc*Area*Skin)/hsc))) * dt   
#dSC=(A-E)*dt
  dSC = ((((Ksc*D_sc*Area*Skin)/hsc)) - ((kelim1nap*SC) + (kelim2nap*SC))) * dt   
#dU=E*dt
  dU = ((kelim1nap*SC) + (kelim2nap*SC)) * dt

  return dSkin, dSC, dU         

############ end of function ############################################## 


for i in range(len(t)-1):   
      
      dSkin, dSC, dU = ODE_dI(i, Ksc, D_sc, Area, Skin[i], hsc, S[i], kelim1nap, kelim2nap, SC[i],dt)  

      Skin[i+1]= Skin[i] + dSkin
      SC[i+1] = SC[i] + dSC
      U[i+1] = U[i] + dU
        
      if (t[i]>180):# stop absorbing after 180 minutes
          Skin[i+1]= 0 

# plot solution
plt.plot(t,Skin, color = "#1f77b4", linewidth = 3) #Skin surface
plt.plot(t,SC,color = "#2ca02c", linewidth = 3) #Stratum corneum
plt.plot(t,U,color = "#ff7f0e", linewidth = 3) #Urine
plt.yscale('log')
plt.xlabel("Time, mins")
plt.ylabel('Naphthalene, ng') 
plt.legend(['Skin surface','Stratum corneum (SC)', 'Urine'])

max_y_skin = max(Skin)
max_x_skin = t[np.argmax(Skin)]
print("Maximum concentration on Skin is", max_y_skin, "ng/g at time", max_x_skin, "minutes")


max_y_SC = max(SC)
max_x_SC = t[np.argmax(SC)]
print("Maximum concentration in Strantum Corneum is", max_y_SC, "ng/g at time", max_x_SC, "minutes")


max_y_U = max(U)
max_x_U = t[np.argmax(U)]
print("Maximum concentration in Urine is", max_y_U,"ng/g at time", max_x_U, "minutes")


plt.savefig("abc.svg", dpi=300, bbox_inches='tight')
files.download("abc.svg")

"""Zoom in of beginning of model, to see the changes in frist 50 minutes of exposure"""

import numpy as np
import matplotlib.pyplot as plt
from google.colab import files


# time period of absorption = 0 - 180 mins (approx time until shower)
# time period of excretion = 43 min fire + 18 hours afterwards = 1123 mins

tmax=1123            #minutes
dt = 20               #time step 
t=np.arange(0,tmax,dt) 

############  Constants ############################################## 
Ksc = 0.177          #partioning coefficient
D_sc= 0.000093/60    #Diffusion coefficient in cm^2/h, , converted to cm^2/min 
Area = 1             #total skin surface area in cm^2
hsc = 0.001485       #skin depth in cm (14.85 um)
kelim1nap= 0.201/60  #urinary elimination coefficient for 1-hydroxynaphthalene in 1/h -> divided by 60=1/min
kelim2nap= 0.288/60  #urinary elimination coefficient for 2-hydroxynaphthalene in 1/h -> divided by 60=1/min
Cskin = 2.36         #max dose in ng/cm2
#########################################################################

S = np.zeros(len(t))
Skin = np.zeros(len(t))
SC = np.zeros(len(t))
U = np.zeros(len(t))

Skin[0]=Cskin  #exposure given in a single dose
SC[0]=0
U[0]=0

############ define function to solve model using Euler #################

def ODE_dI(t, Ksc, D_sc, Area, Skin, hsc, S, kelim1nap, kelim2nap, SC,dt):

#dS=(-A)*dt
  dSkin = - ((((Ksc*D_sc*Area*Skin)/hsc))) * dt   
#dSC=(A-E)*dt
  dSC = ((((Ksc*D_sc*Area*Skin)/hsc)) - ((kelim1nap*SC) + (kelim2nap*SC))) * dt   
#dU=E*dt
  dU = ((kelim1nap*SC) + (kelim2nap*SC)) * dt

  return dSkin, dSC, dU         

############ end of function ############################################## 


for i in range(len(t)-1):   
      
      dSkin, dSC, dU = ODE_dI(i, Ksc, D_sc, Area, Skin[i], hsc, S[i], kelim1nap, kelim2nap, SC[i],dt)  

      Skin[i+1]= Skin[i] + dSkin
      SC[i+1] = SC[i] + dSC
      U[i+1] = U[i] + dU
        
      if (t[i]>180):# stop absorbing after 180 minutes
          Skin[i+1]= 0 

plt.figure(figsize=(10,3))

# plot change in skin concentration from 0-200 mins

plt.subplot(1,2,1)
plt.plot(t,Skin, color = "#1f77b4", linewidth = 3) #Skin surface
plt.xlim(0,200)
plt.ylim(2.27,2.36)
plt.xlabel("Time, mins")
plt.ylabel('Naphthalene, ng') 
plt.legend(['Skin surface'])

# plot change in SC and urine from 0-200 mins

plt.subplot(1,2,2)
plt.plot(t,SC,color = "#2ca02c", linewidth = 3) #Stratum corneum
plt.plot(t,U,color = "#ff7f0e", linewidth = 3) #Urine
plt.xlim(0,50)
plt.ylim(0,.02)
plt.xlabel("Time, mins")
plt.legend(['Stratum \n corneum (SC)', 'Urine'], loc='best', )

plt.tight_layout()

plt.savefig("abc.svg", dpi=300, bbox_inches='tight')
files.download("abc.svg")

"""Investigating the effect of time of decontamination on the absorption of naphthalene"""

import numpy as np
import matplotlib.pyplot as plt
from google.colab import files


# time period of excretion = 43 min fire + 18 hours afterwards = 1123 mins

tmax=1123            #minutes
dt = 20               #time step 
t=np.arange(0,tmax,dt) 

############  Constants ############################################## 
Ksc = 0.177          #partioning coefficient
D_sc= 0.000093/60    #Diffusion coefficient in cm^2/h, converted to cm^2/min 
Area = 1             #skin surface area in cm^2
hsc = 0.001485       #skin depth in cm (14.85 um)
kelim1nap= 0.201/60  #urinary elimination coefficient for 1-hydroxynaphthalene in 1/h, converted to 1/min 
kelim2nap= 0.288/60  #urinary elimination coefficient for 2-hydroxynaphthalene in 1/h, converted to 1/min 
Cskin = 2.36         #max dose in ng/cm2
#########################################################################

S = np.zeros(len(t))
Skin = np.zeros(len(t))
SC = np.zeros(len(t))
U = np.zeros(len(t))

Skin[0]=Cskin  #exposure given in a single dose
SC[0]=0
U[0]=0

############ define function to solve model using Euler #################

def ODE_dI(t, Ksc, D_sc, Area, Skin, hsc, S, kelim1nap, kelim2nap, SC,dt):

#dS=(-A)*dt
  dSkin = - ((((Ksc*D_sc*Area*Skin)/hsc))) * dt   
#dSC=(A-E)*dt
  dSC = ((((Ksc*D_sc*Area*Skin)/hsc)) - ((kelim1nap*SC) + (kelim2nap*SC))) * dt   
#dU=E*dt
  dU = ((kelim1nap*SC) + (kelim2nap*SC)) * dt

  return dSkin, dSC, dU       

############ end of function ############################################## 

showertime=[43,180,600,1123]

plt.figure(figsize=(20,10))

for y in range(len(showertime)):
  for i in range(len(t)-1):
    dSkin, dSC, dU = ODE_dI(i, Ksc, D_sc, Area, Skin[i], hsc, S[i], kelim1nap, kelim2nap, SC[i],dt)  
    Skin[i+1]= Skin[i] + dSkin
    SC[i+1] = SC[i] + dSC
    U[i+1] = U[i] + dU 
    if (t[i]>showertime[y]):
      Skin[i+1]= 0

  max_y_skin = max(Skin)
  max_x_skin = t[np.argmax(Skin)]
  print("Maximum concentration on Skin is", max_y_skin, "ng/g at time", max_x_skin, "minutes, after waiting", showertime[y], "mins before decon")

  max_y_SC = max(SC)
  max_x_SC = t[np.argmax(SC)]
  print("Maximum concentration in Strantum Corneum is", max_y_SC, "ng/g at time", max_x_SC, "minutes, after waiting", showertime[y], "mins before decon" )

  max_y_U = max(U)
  max_x_U = t[np.argmax(U)]
  print("Maximum concentration in Urine is", max_y_U,"ng/g at time", max_x_U, "minutes, after waiting", showertime[y], "mins before decon")


  plt.subplot(3,2,y+1)
  plt.plot(t,Skin, color = "#1f77b4", linewidth = 3) #Skin surface
  plt.plot(t,SC,color = "#2ca02c", linewidth = 3) #Stratum corneum
  plt.plot(t,U,color = "#ff7f0e", linewidth = 3) #Urine
  plt.xlabel("Time, mins", fontsize = 15)
  plt.ylabel('Naphthalene, ng', fontsize = 15)
  plt.yscale('log') 
  plt.ylim(0,3)
  plt.title("Time of decon = "+str(showertime[y]) +" mins", fontsize = 20)


plt.figlegend(['Skin surface','Stratum corneum (SC)', 'Urine'], loc="center",bbox_to_anchor=(0.25, -0.12,0.5, 0.5), ncol = 3, fontsize = 25) #x ,y , width, length
plt.suptitle('Naphthalene absorption with respect to time of decon', fontsize=30, x=0.5, y=1.04)
plt.tight_layout()

print(U)

print(t)

plt.savefig("abc.svg", dpi=300, bbox_inches='tight')
files.download("abc.svg")

"""Investigating differences in the absorption of naphthalene at different locations on the body"""

import numpy as np
import matplotlib.pyplot as plt
from google.colab import files

# time period of excretion = 43 min fire + 18 hours afterwards = 1123 mins

tmax=1123            #minutes
dt = 20               #time step 
t=np.arange(0,tmax,dt) 

############  Constants ############################################## 
Ksc = 0.177           #partioning coefficient
D_sc= 0.000093/60     #Diffusion coefficient in cm^2/h, converted to cm/min
Area = 1              #total skin surface area in cm^2
kelim1nap= 0.201/60   #urinary elimination coefficient for 1-hydroxynaphthalene in 1/h, converted to 1/min 
kelim2nap= 0.288/60   #urinary elimination coefficient for 2-hydroxynaphthalene in 1/h, converted to 1/min 
Cskin = 2.36          #max dose in ng/cm2
#########################################################################

S = np.zeros(len(t))
Skin = np.zeros(len(t))
SC = np.zeros(len(t))
U = np.zeros(len(t))

Skin[0]=Cskin  #exposure given in a single dose
SC[0]=0
U[0]=0

############ define function to solve model using Euler #################

def ODE_dI(t, Ksc, D_sc, Area, Skin, hsc, S, kelim1nap, kelim2nap, SC,dt):

#dS=(-A)*dt
  dSkin = - ((((Ksc*D_sc*Area*Skin)/hsc))) * dt   
#dSC=(A-E)*dt
  dSC = ((((Ksc*D_sc*Area*Skin)/hsc)) - ((kelim1nap*SC) + (kelim2nap*SC))) * dt   
#dU=E*dt
  dU = ((kelim1nap*SC) + (kelim2nap*SC)) * dt

  return dSkin, dSC, dU       

############ end of function ############################################## 

plt.figure(figsize=(20,10))

hsc = [0.00093, 0.00109, 0.00062, 0.00063] #skin depth in cm from https://www.medicaljournals.se/acta/download/10.2340/00015555-0875/
body_position = ["Back of Hand", "Outer Forearm", "Inner Forearm", "Temple"]


for y in range(len(hsc)) and range(len(body_position)):
  for i in range(len(t)-1):
    dSkin, dSC, dU = ODE_dI(i, Ksc, D_sc, Area, Skin[i], hsc[y], S[i], kelim1nap, kelim2nap, SC[i],dt)  
    Skin[i+1]= Skin[i] + dSkin
    SC[i+1] = SC[i] + dSC
    U[i+1] = U[i] + dU
    if (t[i]>180):
      Skin[i+1]= 0
  
  max_y_skin = max(Skin)
  max_x_skin = t[np.argmax(Skin)]
  print("Maximum concentration on Skin is", max_y_skin, "ng/g at time", max_x_skin, "minutes, on", body_position[y])

  max_y_SC = max(SC)
  max_x_SC = t[np.argmax(SC)]
  print("Maximum concentration in Strantum Corneum is", max_y_SC, "ng/g at time", max_x_SC, "minutes, on", body_position[y])


  max_y_U = max(U)
  max_x_U = t[np.argmax(U)]
  print("Maximum concentration in Urine is", max_y_U,"ng/g at time", max_x_U, "minutes, on", body_position[y])
 

  plt.subplot(3,2,y+1)
  plt.plot(t,Skin, color = "#1f77b4", linewidth = 3) #Skin surface
  plt.plot(t,SC,color = "#2ca02c", linewidth = 3) #Stratum corneum
  plt.plot(t,U,color = "#ff7f0e", linewidth = 3) #Urine
  plt.xlabel("Time, mins", fontsize = 15)
  plt.ylabel('Naphthalene, ng', fontsize = 15) 
  plt.yscale('log')
  plt.ylim(0,4)
  plt.title(str(body_position[y]), fontsize = 20)

plt.figlegend(['Skin surface','Stratum corneum (SC)', 'Urine'], loc="center",bbox_to_anchor=(0.25, -0.12,0.5, 0.5), ncol = 3, fontsize = 25) #x ,y , width, length
plt.suptitle('Naphthalene absorption by body location', fontsize=30, x=0.5, y=1.04)
plt.tight_layout()

plt.savefig("abc.svg", dpi=300, bbox_inches='tight')
files.download("abc.svg")

"""Investigating differences in the absorption of naphthalene at different locations on the body and time of decon"""

import numpy as np
import matplotlib.pyplot as plt
import itertools
from google.colab import files

# time period of excretion = 43 min fire + 18 hours afterwards = 1123 mins

tmax=1123            #minutes
dt = 20               #time step 
t=np.arange(0,tmax,dt) 

############  Constants ############################################## 
Ksc = 0.177           #partioning coefficient
D_sc= 0.000093/60     #Diffusion coefficient in cm^2/h, , converted to cm/min
Area = 1              #total skin surface area in cm^2
kelim1nap= 0.201/60   #urinary elimination coefficient for 1-hydroxynaphthalene in 1/h -> divided by 60=1/min
kelim2nap= 0.288/60   #urinary elimination coefficient for 2-hydroxynaphthalene in 1/h -> divided by 60=1/min
Cskin = 2.36          #max dose in ng/cm2
#########################################################################

S = np.zeros(len(t))
Skin = np.zeros(len(t))
SC = np.zeros(len(t))
U = np.zeros(len(t))

Skin[0]=Cskin  #exposure given in a single dose
SC[0]=0
U[0]=0

############ define function to solve model using Euler #################

def ODE_dI(t, Ksc, D_sc, Area, Skin, hsc, S, kelim1nap, kelim2nap, SC,dt):

#dS=(-A)*dt
  dSkin = - ((((Ksc*D_sc*Area*Skin)/hsc))) * dt   
#dSC=(A-E)*dt
  dSC = ((((Ksc*D_sc*Area*Skin)/hsc)) - ((kelim1nap*SC) + (kelim2nap*SC))) * dt   
#dU=E*dt
  dU = ((kelim1nap*SC) + (kelim2nap*SC)) * dt

  return dSkin, dSC, dU       

############ end of function ############################################## 

showertime=[43,180,600,1123]

hsc = [0.00093, 0.00109, 0.00062, 0.00063] #skin depth in cm from https://www.medicaljournals.se/acta/download/10.2340/00015555-0875/
body_position = ["Back of Hand", "Outer Forearm", "Inner Forearm", "Temple"]


combinations = list(itertools.product(showertime, hsc)) #make a list of list with the different combinations of temp+thickness


all = list.append(combinations,(1123, 0.00093))


plt.figure(figsize=(15, 10))

plot_position = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]


body = ["Back of Hand", "Outer Forearm", "Inner Forearm", "Temple","Back of Hand", "Outer Forearm", "Inner Forearm", "Temple","Back of Hand", "Outer Forearm", "Inner Forearm", "Temple","Back of Hand", "Outer Forearm", "Inner Forearm", "Temple"]
for j,k,l,m in zip(combinations,combinations[1:],plot_position,body):
  temp = j[0]
  for i in range(len(t)-1):
    dSkin, dSC, dU = ODE_dI(i, Ksc, D_sc, Area, Skin[i], k[1], S[i], kelim1nap, kelim2nap, SC[i],dt)  
    Skin[i+1]= Skin[i] + dSkin
    SC[i+1] = SC[i] + dSC
    U[i+1] = U[i] + dU 
    if (t[i]>temp):
      Skin[i+1]= 0

  max_y_skin = max(Skin)
  max_x_skin = t[np.argmax(Skin)]
  print("Maximum concentration on Skin is", max_y_skin, "ng/g at time", max_x_skin, "minutes, on", m, "when you wait " + str(temp) + " min to decon")

  max_y_SC = max(SC)
  max_x_SC = t[np.argmax(SC)]
  print("Maximum concentration in Strantum Corneum is", max_y_SC, "ng/g at time", max_x_SC, "minutes, on", m, "when you wait " + str(temp) + " min to decon")


  max_y_U = max(U)
  max_x_U = t[np.argmax(U)]
  print("Maximum concentration in Urine is", max_y_U,"ng/g at time", max_x_U, "minutes, on", m, "when you wait " + str(temp) + " min to decon")


  plt.subplot(5,4,l)
  plt.plot(t,Skin, color = "#1f77b4", linewidth = 3) #Skin surface
  plt.plot(t,SC,color = "#2ca02c", linewidth = 3) #Stratum corneum
  plt.plot(t,U,color = "#ff7f0e", linewidth = 3) #Urine
  plt.xlabel("Time, mins")
  plt.ylabel('Naphthalene, ng') 
 # plt.ylim(0,4)
  plt.yscale('log')
  plt.title(m + "\n" + "Time of decon = "+ str(temp) +" mins")

plt.figlegend(['Skin surface','Stratum corneum (SC)', 'Urine'], loc="center",bbox_to_anchor=(0.25, -0.17,0.5, 0.5), ncol = 3, fontsize = 25) #x ,y , width, length
plt.suptitle('Naphthalene absorption with respect to decon time and body location', fontsize=30, x=0.5, y=1.05)
plt.tight_layout()

plt.savefig("abc.svg", dpi=300, bbox_inches='tight')
files.download("abc.svg")
