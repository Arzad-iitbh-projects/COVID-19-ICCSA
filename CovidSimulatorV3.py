import random
from matplotlib import pyplot as plt
import scipy.integrate as integrate
import numpy as np
import math

Num_Pub_Fac=20     #Number of Public Facilities
Num_Pri_Fac=1
Num_Qua_Fac=1
Fac={}            #list of facilities
service=[]        #probability of person leaving a facility 
routingProbability={}  #probability of person going to each facility given
                       #the facility that he/she is currently in
Pstates={}
Stimestamp={}
Itimestamp={}
rates=[]
lengt=30

slot_length = 1 #minute
average_infection_time = 15 #minutes
average_detection_time = 50*24*60 #minutes
average_recovery_time = 14*24*60 #minutes
average_service_time_pub = 15 #minutes
average_service_time_pri = 5*60 #minutes
average_quarantine_time = 14*24*60 #minutes


simulation_time = 30*24*60 #minutes

Percent_Initial_Infectious= 1 #percentage. it is not fraction.
Pspread=slot_length/average_infection_time
Pdetected=slot_length/average_detection_time
Precovery=slot_length/average_recovery_time

mu_Priv = slot_length/average_service_time_pri
mu_Public = slot_length/average_service_time_pub
theta = 0.8
Pquar = slot_length/average_quarantine_time
t_sim=int(simulation_time/slot_length)
t_settle=int(0.2*simulation_time/slot_length) #to allow the initial dynamics to settle down. no infectious
                                                #individuals in this duration. this ensures a steady state
                                                #is achieved so that the arrival rates into different facilities
                                                #have converged

N=[]    #list of individuals
Inf=[]
Det=[]
N_Temp=[]

def Initialize_Logging():
    for i in range(Num_Pub_Fac+Num_Pri_Fac+Num_Qua_Fac):  #here we update values into Fac,service,RoutingProbability
        Stimestamp.update({i:[]})
        Itimestamp.update({i:[]})
    
def Initialize_Model():
    for i in range(Num_Pub_Fac):  #here we update values into Fac,service,RoutingProbability
        Fac.update({i:[]})
        service.append(mu_Public)
        routingProbability.update({i:[]})
        for j in range(Num_Pub_Fac):
            routingProbability[i].append((1-theta)/(Num_Pub_Fac-1))
        for j in range(Num_Pri_Fac):
            routingProbability[i].append(theta/Num_Pri_Fac)
        routingProbability[i][i]=0.0

    for i in range(Num_Pub_Fac,Num_Pub_Fac+Num_Pri_Fac,1):  #here we update values into Fac,service,RoutingProbability
        Fac.update({i:[]})
        service.append(mu_Priv)
        routingProbability.update({i:[]})
        for j in range(Num_Pub_Fac):
            routingProbability[i].append(1/(Num_Pub_Fac))
        for j in range(Num_Pri_Fac):
            routingProbability[i].append(0.0)
        routingProbability[i][i]=0.0

    for i in range(Num_Pub_Fac+Num_Pri_Fac,Num_Pub_Fac+Num_Pri_Fac+Num_Qua_Fac,1):  #here we update values into Fac,service,RoutingProbability
        Fac.update({i:[]})
        service.append(Pquar)

def Initialize_Population():
    for i in range(Population):  #initial state
        N.append(i)
        N_Temp.append(i)
        Pstates.update({i:['S',i%(Num_Pub_Fac+Num_Pri_Fac)]})
        Fac[Pstates[i][1]].append(i)

def Initialize_Infectious():
    for i in range(math.floor(Population*Percent_Initial_Infectious/100)):
        p=random.choices(N_Temp,k=1)
        Inf.append(p[0])
        N_Temp.remove(p[0])
        Pstates[p[0]][0]='I'

def position(name):   #this function tells us which facility a given person is in
    return(Pstates[name][1])

def rotation(daynum,update_flag): #this function tells us the movement of people between facilities
    for i in range(Population):
        j=position(i)
        p=random.random()
        if p<=service[j]:            #probability that person is leaving the facility
            if j<(Num_Pub_Fac+Num_Pri_Fac):
                Fac[j].remove(i)
                target=random.choices(range(Num_Pub_Fac+Num_Pri_Fac),routingProbability[j])
                Pstates[i][1]=target[0]
                Fac[Pstates[i][1]].append(i)
                if update_flag:
                    if Pstates[i][0]=='S':
                        Stimestamp[Pstates[i][1]].append(daynum)
                    if Pstates[i][0]=='I':
                        Itimestamp[Pstates[i][1]].append(daynum)
                
def infectionSpread():
    for key in Fac:           #goes through each facility
        infected=False
        for i in Fac[key]:    #goes through people in each facility
            if infected==False and Pstates[i][1]<Num_Pub_Fac:
                if Pstates[i][0]=='I':
                    infected=True
                    break
        if infected==True:
            for j in Fac[key]:
                if Pstates[j][0]=='S':
                    p=random.random()
                    if p<=Pspread:
                        Pstates[j][0]='I'
                        Inf.append(j)

def detected():
    for i in Inf:
        p=random.random()
        if p<=Pdetected:
            Pstates[i][0]='D'
            Inf.remove(i)
            Det.append(i)
            k=position(i)
            Fac[k].remove(i)
            Fac[Num_Pub_Fac+Num_Pri_Fac].append(i)
            Pstates[i][1]=Num_Pub_Fac+Num_Pri_Fac

def recovery():
    for i in Inf:
        p=random.random()
        if p<=Precovery:
            Inf.remove(i)
            Pstates[i][0]='R'

Initialize_Model()
Population=int(input('Enter population: '))

Initialize_Population()
#all the individuals are susceptible now. no infection.
#this ensures that the initial transient for the arrival rates
#in the system has died out and we have converged to a steady state.
for k in range(t_settle):
    rotation(k,0) #call rotation with logging disabled

g=[]
Initialize_Infectious()
Initialize_Logging()    #Initialize the data structures for logging
for k in range(t_sim):
    infectionSpread()
    detected()
    recovery()
    rotation(k,1) #call rotation with logging enabled
    g.append(len(Inf))
    
figure,axis = plt.subplots(1,2)
axis[0].set_xlabel('Time (in minutes)')
axis[1].set_xlabel('Time (in minutes)')
axis[0].set_ylabel('Susceptible Arrival Rates')
axis[1].set_ylabel('Infected Arrival Rates')

for j in range(Num_Pub_Fac+Num_Pri_Fac):
    Srates=[]
    labelled=False
    for i in range(lengt,len(Stimestamp[j]),1):
        Srates.append(min(i,lengt)/(slot_length*(Stimestamp[j][i]-Stimestamp[j][max(0,i-lengt)])))    
    if j==0:
        axis[0].plot(Stimestamp[j][lengt:],Srates,label='Public Facility #1')
    if j==Num_Pub_Fac+Num_Pri_Fac-1:
        axis[0].plot(Stimestamp[j][lengt:],Srates,label='Private Facility')
        
       
    


for j in range(Num_Pub_Fac+Num_Pri_Fac):
    rates=[]
    for i in range(lengt,len(Itimestamp[j]),1):
        rates.append(min(i,lengt)/(slot_length*(Itimestamp[j][i]-Itimestamp[j][max(0,i-lengt)])))
    if j==0:
        axis[1].plot(Itimestamp[j][lengt:],rates,label='Public Facility #1')
    if j==Num_Pub_Fac+Num_Pri_Fac-1:
        axis[1].plot(Itimestamp[j][lengt:],rates,label='Private Facility')
        


axis[0].legend()
axis[1].legend()
plt.show()

