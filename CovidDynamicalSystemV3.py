import random
from matplotlib import pyplot as plt
import scipy.integrate as integrate
import numpy as np
Num_Pub_Fac=50
Num_Pri_Fac=1

def p_s(lamb,mu,rho):

    #lamb        =1/1000   #denominator is the mean interarrival time in minutes
    #mu          =1/100   #denominator is the mean service requirement, in minutes
    #rho         =1/300   #denominator is the mean time it takes to infect someone, in minutes

    result = integrate.quad(lambda t: np.exp(
                                    -(rho+mu)*t
                                    - (lamb/mu)*(1- np.exp(-mu*t))
                                    ), 0, 600)
    
    B_hat_ = 1 + (1/lamb)*(rho+mu -1/result[0])
    ExpB    = (np.exp(lamb/mu)-1)/lamb
    #print(lamb, mu, ExpB)
    # print(ExpB)
    # p=0 = C + D p=1
    # p=1 = E p=0 + F
    # p=0 = (C+DF)/(1-DE)
    # p>=1 = A p=0 + B

    A = (1-B_hat_)/((rho+mu)*ExpB)
    B = (mu/(rho+mu))*( 1 - (1-B_hat_)/((rho+mu)*ExpB) )
    C = mu/(lamb+mu)
    D = lamb/(mu+lamb)
    E = B_hat_
    F = (1-B_hat_)*mu/((rho+mu))

    p_eq_0 = (C+D*F)/(1-D*E)
    p_geq_1 = A*p_eq_0+B

    pi_0 = np.exp(-lamb/mu)

    ps = pi_0*p_eq_0 + (1-pi_0)*p_geq_1
        
    return(ps)

lambdas=[]#arrival rates of the infectious individuals
omegas=[]#arrival rates of the susceptible individuals
mus=[]
rhos=[]
lambdas_temp=[]
omegas_temp=[]

lamb_evo=[]
omega_evo=[]
routingProbability={}  #probability of person going to each facility given
                       #the facility that he/she is currently in

average_infection_time = 15 #minutes
average_detection_time = 50*24*60 #minutes
average_recovery_time = 14*24*60 #minutes
average_service_time_pub = 15 #minutes
average_service_time_pri = 5*60 #minutes
average_quarantine_time = 14*24*60 #minutes
t=100 #number of iterations

percentage_initial_infectious=1 #percentage. not fraction.

zeta = 1/average_detection_time + 1/average_recovery_time
rho_Public=1/average_infection_time
rho_Priv=0.0
mu_Priv = 1/average_service_time_pri
mu_Public = 1/average_service_time_pub
theta = 0.8
delta=[]


for i in range(Num_Pub_Fac):  #here we update values into Fac,service,RoutingProbability
    delta.append(zeta/(zeta+1/(average_service_time_pub))) #the probability of infectious node becoming recovered or detected while in the facility
    routingProbability.update({i:[]})
    for j in range(Num_Pub_Fac):
        routingProbability[i].append((1-theta)/(Num_Pub_Fac-1))
    for j in range(Num_Pri_Fac):
        routingProbability[i].append(theta/Num_Pri_Fac)
    routingProbability[i][i]=0.0

for i in range(Num_Pub_Fac,Num_Pub_Fac+Num_Pri_Fac,1):  #here we update values into Fac,service,RoutingProbability
    delta.append(zeta/(zeta+1/(average_service_time_pri)))#the probability of infectious node becoming recovered or detected while in the facility
    routingProbability.update({i:[]})
    for j in range(Num_Pub_Fac):
        routingProbability[i].append(1/(Num_Pub_Fac))
    for j in range(Num_Pri_Fac):
        routingProbability[i].append(0.0)
    routingProbability[i][i]=0.0

Population=int(input('Enter population: '))
#the following part computes the relation between the number of
#individuals in the system and the corresponding overall arrival rates.
#this assumes that the individuals are all at the same private facility initially.
for i in range(Num_Pub_Fac):
    omegas.append(mu_Priv*Population/(Num_Pub_Fac))
    lambdas.append(0.2/10)
    lambdas_temp.append(0.0)
    omegas_temp.append(0.0)
    mus.append(mu_Public)
    rhos.append(rho_Public)

for i in range(Num_Pub_Fac,Num_Pub_Fac+Num_Pri_Fac,1):
    omegas.append(0.00001)
    lambdas.append(4/10)
    lambdas_temp.append(0.0)
    omegas_temp.append(0.0)
    mus.append(mu_Priv)
    rhos.append(rho_Priv)


for time in range(int(0.2*t)):
    for i in range(Num_Pub_Fac+Num_Pri_Fac):
        omegas_temp[i]=0.0
        for j in range(Num_Pub_Fac+Num_Pri_Fac):
            if j!=i:
                omegas_temp[i]=omegas_temp[i]+(omegas[j]*routingProbability[j][i])
    for i in range(Num_Pub_Fac+Num_Pri_Fac):
        omegas[i]=omegas_temp[i]

for i in range(Num_Pub_Fac+Num_Pri_Fac):
    lambdas[i]=(percentage_initial_infectious/100)*omegas[i]
    omegas[i]=(1-percentage_initial_infectious/100)*omegas[i]

for time in range(t):
    lamb_evo.append([])
    omega_evo.append([])
    for i in range(Num_Pub_Fac+Num_Pri_Fac):
        lamb_evo[-1].append(lambdas[i])
        omega_evo[-1].append(omegas[i])
        lambdas_temp[i] = 0.0
        omegas_temp[i]=0.0
        for j in range(Num_Pub_Fac+Num_Pri_Fac):
            if j!=i:
                omegas_temp[i]=omegas_temp[i]+(omegas[j]*p_s(lambdas[j],mus[j],rhos[j])*routingProbability[j][i])
                lambdas_temp[i]=lambdas_temp[i]+(omegas[j]*(1-p_s(lambdas[j],mus[j],rhos[j]))*routingProbability[j][i] + lambdas[j]*(1-delta[j])*routingProbability[j][i])

    for i in range(Num_Pub_Fac+Num_Pri_Fac):
        lambdas[i]=lambdas_temp[i]
        omegas[i]=omegas_temp[i]

fig,axis=plt.subplots(1,2)
axis[0].set_xlabel('Iteration')
axis[1].set_xlabel('Iteration')
axis[0].set_ylabel('Susceptible Arrival Rates')
axis[1].set_ylabel('Infected Arrival Rates')

for j in range(len(omega_evo[0])):
    k=[]
    for i in range(len(omega_evo)):
        k.append(omega_evo[i][j])
    if j==0:
        axis[0].plot(k,label='Public Facility #1')
    if j==Num_Pub_Fac+Num_Pri_Fac-1:
        axis[0].plot(k,label='Private Facility')

for j in range(len(lamb_evo[0])):
    k=[]
    for i in range(len(lamb_evo)):
        k.append(lamb_evo[i][j])
    if j==0:
        axis[1].plot(k,label='Public Facility #1')
    if j==Num_Pub_Fac+Num_Pri_Fac-1:
        axis[1].plot(k,label='Private Facility')

axis[0].legend()
axis[1].legend()
plt.show()





