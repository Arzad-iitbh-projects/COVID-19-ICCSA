import csv
import random
import scipy.integrate as integrate
import numpy as np
from matplotlib import pyplot as plt

Ts=1          #one time slot is of length 1ms
a_rate=Ts/1000    #denominator is the mean interarrival time in seconds
mu=Ts/100       #denominator is the mean service requirement in seconds
p_spread=Ts/300  #denominator is the mean time (in seconds) that it takes for
                #a susceptible individual to become infected when in contact


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

    # p=0 = C + D p=1
    # p=1 = E p=0 + F
    # p=0 = (C+DF)/(1-DE)
    # p>=1 = A p=0 + B

    #A = (rho/(mu+rho))*(1-B_hat_)/((rho+mu)*ExpB)
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

def get_estimate(a_rate, mu, p_spread):
    susc=[]
    inf=[]
    num_susceptible_unchanged=0
    num_susceptible_arrivals=0

    not_empty=0
    not_empty_direct=0
    num_slots=0

    while num_susceptible_arrivals<20000:
        if len(inf):
            not_empty_direct=not_empty_direct+1
        num_slots=num_slots+1
        #t in range(10000000):
        #arrival
        if random.random()<a_rate:
            inf.append('I')
            #print(" arrival of inf ",inf)
        if random.random()<a_rate:
    #        if len(inf):
            susc.append('S')
            #print(" arrival of susc ",susc)
            num_susceptible_arrivals=num_susceptible_arrivals+1
            if len(inf):
                not_empty=not_empty+1;

        #infection
        if len(inf):
            for i in susc:
                if random.random()<p_spread:
                    #print(" infection spread initiation ",susc,inf)
                    susc.remove('S')
                    #inf.append('I')    #so that the changed infectious people do not contribute to the spread
                    #print(" infection spread completion ",susc,inf)

        #departure
        for i in susc:
            if random.random()<mu:
                susc.remove('S')
                num_susceptible_unchanged=num_susceptible_unchanged+1
        for i in inf:
            if random.random()<mu:
                inf.remove('I')

    return(num_susceptible_unchanged/num_susceptible_arrivals)

#csvfile=open('result.csv', 'w', newline='')
#data=[("lambda","mu","rho","p_s sim","p_s num")]
#for row in data:
#    line = ','.join(row)
#    csvfile.write(line + '\n')
#csvfile.close()

fig,axis=plt.subplots()
axis.set_xlabel('p_spread')
axis.set_ylabel('p_s')

result={}
result["p_spread"]=[]
result["lambda"]=[]
result["mu"]=[]
result["p_s_sim"]=[]
result["p_s_num"]=[]


for p_spread_inv in range(10,1000,20):
    for lambdainv in range(100,1000,1200):
        for muinv in range(10,100,120):
            p_s_sim=get_estimate(1.0/lambdainv, 1.0/muinv, 1.0/p_spread_inv)
            p_s_num=p_s(1.0/lambdainv,1.0/muinv,1.0/p_spread_inv)
            result["p_spread"].append(1.0/p_spread_inv)
            result["lambda"].append(1.0/lambdainv)
            result["mu"].append(1.0/muinv)
            result["p_s_sim"].append(p_s_sim)
            result["p_s_num"].append(p_s_num)
            print("{0:.4f}".format(1.0/lambdainv),"{0:.4f}".format(1.0/muinv),"{0:.4f}".format(1.0/p_spread_inv),"{0:.4f}".format(p_s_sim),"{0:.4f}".format(p_s_num))
            #csvfile=open('result.csv', 'a', newline='')
            #data=[(str("{0:.4f}".format(1.0/lambdainv)),str("{0:.4f}".format(1.0/muinv)),str("{0:.4f}".format(1.0/p_spread_inv)),str("{0:.4f}".format(p_s_sim)),str("{0:.4f}".format(p_s_num)))]
            #for row in data:
            #    line = ','.join(row)
            #    csvfile.write(line + '\n')
            #csvfile.close()
#    axis.scatter(1.0/p_spread_inv,p_s_sim,label='p_s_sim')
#    axis.scatter(1.0/p_spread_inv,p_s_num,label='p_s_num')
            #mywriter.writerow()
axis.legend()

axis.plot(result["p_spread"],result["p_s_sim"],label="p_s_sim")
axis.plot(result["p_spread"],result["p_s_num"],label="p_s_num")
plt.show()

#    print(" Probability of system not empty (simulation) ", "{0:.2f}".format(not_empty/num_susceptible_arrivals))
#    print(" Probability of system not empty (analytical) ", "{0:.2f}".format(1-numpy.exp(-a_rate/mu)))
#    print(" Probability of system not empty (simu... direct) ", "{0:.2f}".format(not_empty_direct/num_slots))
