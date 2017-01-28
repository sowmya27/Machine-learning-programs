__author__ = 'annapurnaannadatha'
import numpy as np


def estimate_T(P):
    T1 = P[0].sum()/5
    T2 = 1 - T1
    T=[T1,T2]
    return T

def estimate_O(P,X):
    U=[]
    for j in range(0,len(O)):
        U1 = 0
        for i in range(0,len(X)):
            U1 = U1 + P[j][i]*X[i]
        U1 = U1/(P[j].sum()*10)
        U.append(U1)
    return U


def estep(O,T,X):
    P = np.empty(shape = (2,5))
    for j in range(0,len(O)):
        for i in range(0,len(X)):
            f = T[j] * O[j]**X[i] *(1-O[j])**(10-X[i])/(T[0]* O[0]**X[i] *(1-O[0])**(10-X[i]) + T[1]* O[1]**X[i] *(1-O[1])**(10-X[i]) )
            P[j][i]= f
    return P


O = [0.3,0.9]
T = [0.2,0.8]
X = [8,5,9,4,7]
Iterations = int(input("enter no of iterations:"))
print("Intial O=", O)
print("Intial T=",T)
print("Intial X=",X)

for i in range(0,Iterations):
    P = estep(O,T,X)
    T = estimate_T(P)
    O = estimate_O(P,X)
    print("Iteration no:",i+1)
    print("O=", O)
    print("T=",T)
    print()

print()
print()


'''
T = [0.7593,0.2407]
O = [0.6918,0.5597]
P = estep(O,T,X)
print(P)
'''




