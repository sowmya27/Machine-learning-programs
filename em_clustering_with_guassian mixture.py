__author__ = 'annapurnaannadatha'
#em clustering with guassian mixture
import numpy as np

def estep(X,U1,U2,S1,S2,T):
    P = np.empty(shape = (2,20))
    D1 = S1[0][1]/np.sqrt(S1[0][0]*S1[1][1])
    D2 = S2[0][1]/np.sqrt(S2[0][0]*S2[1][1])
    print("D1=",D1)
    print("D2=",D2)
    for i in range(0,20):
        z1 = ((X[i][0]-U1[0][0])**2/(S1[0][0]) - 2*D1*(X[i][0]-U1[0][0])*(X[i][1]-U1[1][0])/(S1[0][1]/D1) + (X[i][1]-U1[1][0])**2/S1[1][1])/(1-D1**2)
        f1 = np.exp(-z1/2)/(2*(np.pi)*(S1[0][1]/D1)*np.sqrt(1-D1**2))

        z2 = ((X[i][0]-U2[0][0])**2/(S2[0][0]) - 2*D2*(X[i][0]-U2[0][0])*(X[i][1]-U2[1][0])/(S2[0][1]/D2) + (X[i][1]-U2[1][0])**2/S2[1][1])/(1-D2**2)
        f2 = np.exp(-z2/2)/(2*(np.pi)*(S2[0][1]/D2)* np.sqrt(1-D2**2))

        P[0][i] = T[0][0] * f1 /(T[0][0]*f1 + T[0][1]*f2)
        print("P[1][%d]="%(i+1),P[0][i])

        P[1][i] = T[0][1] * f2 /(T[0][0]*f1 + T[0][1]*f2)
        print("P[2][%d]="%(i+1),P[1][i])
    return P

def estimate_T(P):
    T1 = P[0].sum()/20
    T2 = 1 - T1
    T = np.array([[T1,T2]])
    return T

def estimate_U1(P,X):
    Ux =0
    Uy =0
    for i in range(0,20):
        Ux = Ux + P[0][i]*X[i][0]
    Ux = Ux/(P[0].sum())
    for i in range(0,20):
        Uy = Uy + P[0][i]*X[i][1]
    Uy = Uy/(P[0].sum())
    U = np.array([[Ux],[Uy]])
    return U

def estimate_U2(P,X):
    Ux =0
    Uy =0
    for i in range(0,20):
        Ux = Ux + P[1][i]*X[i][0]
    Ux = Ux/(P[1].sum())
    for i in range(0,20):
        Uy = Uy + P[1][i]*X[i][1]
    Uy = Uy/(P[1].sum())
    U = np.array([[Ux],[Uy]])
    return U

def estimate_S1(P,X,U1):
    sx =np.zeros(shape=(2,2))
    for i in range(0,20):
        xu = np.array([[X[i][0]-U1[0][0]],[X[i][1]-U1[1][0]]])
        sx = sx + P[0][i] * xu * xu.T
    S = sx/(P[0].sum())
    return S

def estimate_S2(P,X,U2):
    sx =np.zeros(shape=(2,2))
    for i in range(0,20):
        xu = np.array([[X[i][0]-U2[0][0]],[X[i][1]-U2[1][0]]])
        sx = sx + P[1][i] * xu * xu.T
    S = sx/(P[1].sum())
    return S

X = np.array([[3.60,79],
              [1.800,54],
              [2.283,62],
              [3.333,74],
              [2.883,55],
              [4.533,85],
              [1.950,51],
              [1.833,54],
              [4.700,88],
              [3.600,85],
              [1.600,52],
              [4.350,85],
              [3.917,84],
              [4.200,78],
              [1.750,62],
              [1.800,51],
              [4.700,83],
              [2.167,52],
              [4.800,84],
              [1.750,47]
             ])




U1 = np.array([[2.5],
               [65.0]])
U2 = np.array([[3.5],
               [70.0]])
S1 = np.array([[1.0,5.0],
               [5.0,100.0]])
S2 = np.array([[2.0,10.0],
               [10.0,200.0]])
T = np.array([[0.6,0.4]])

for i in range(0,100):
    print("Iteration number :",i+1)
    P = estep(X,U1,U2,S1,S2,T)
    #print(P)
    print()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("After Iteration",i+1)
    print("Restimating values....")

    T = estimate_T(P)
    print("T=")
    print(T)
    print()

    U1 = estimate_U1(P,X)
    print("U1=")
    print(U1)
    print()

    U2 = estimate_U2(P,X)
    print("U2=")
    print(U2)
    print()

    S1 = estimate_S1(P,X,U1)
    print("S1=")
    print(S1)
    print()

    S2 = estimate_S2(P,X,U2)
    print("S2=")
    print(S2)
    print()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
