# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:33:10 2018

@author: USUARIO
"""
#%% Librerías que se utilizan

import numpy as np
from sympy import Matrix,sin,cos
from scipy.linalg import solve
import matplotlib.pyplot as plt
import meshio  
from mpl_toolkits.mplot3d import Axes3D


#%% Se importa el archivo .msh que se va a analizar y se definen parámateros de
#   la ecuación.
points, cells, point_data, cell_data, field_data = meshio.read('cuadrado.msh')
g = cells['line']
condi = list(set(g.flatten()))
w = 1
v = 1
Ka = w**2/v**2

Pos = points
nodos = int (np.size(Pos)/3)

#%% Definición de la función de cargas.

def f(x,y):
    f = x**2+y**2
    return f


#%% Definición de la función que realiza la cuadratura gaussiana con 4 puntos.

def Intgauss(e,n,W1,W2):
    H = Matrix([[1/4*(1-e)*(1-n), 1/4*(1+e)*(1-n), 1/4*(1+e)*(1+n), 1/4*(1-e)*(1+n)]])
    Ben = Matrix([[0.25*n - 0.25, -0.25*n + 0.25, 0.25*n + 0.25, -0.25*n - 0.25],\
          [0.25*e - 0.25, -0.25*e - 0.25, 0.25*e + 0.25, -0.25*e + 0.25]])
    Pos3 = np.zeros((4,2))
    for j in range (4):
        a = int(Con2[0][j])
        x = Pos[a,0]
        y = Pos[a,1]
        Pos3[j][0]= x
        Pos3[j][1]= y
        
    Jaux = Ben*Matrix([[Pos3[0][0],Pos3[0][1]],[Pos3[1][0],Pos3[1][1]],[Pos3[2][0],\
                    Pos3[2][1]],[Pos3[3][0],Pos3[3][1]]])
    J = np.zeros((2,2))
    for i in range (2):
        J[0][i]=Jaux[0,i]
        J[1][i]=Jaux[1,i]
    Jdet = np.linalg.det(J)
    Jinv = np.linalg.inv(J)
    B = Jinv*Ben
    Bt = B.transpose()
    Ht = H.transpose()
    a = Bt*B*Jdet*W1*W2
    b = Ht*H*Jdet*W1*W2
    return a,b
        

#%% Función en donde se encuentran los vectores de carga elementales.

def gaussbe(e,n,W1,W2):
    Ht = [1/4*(1-e)*(1-n), 1/4*(1+e)*(1-n), 1/4*(1+e)*(1+n), 1/4*(1-e)*(1+n)]
    Ben = Matrix([[0.25*n - 0.25, -0.25*n + 0.25, 0.25*n + 0.25, -0.25*n - 0.25],\
          [0.25*e - 0.25, -0.25*e - 0.25, 0.25*e + 0.25, -0.25*e + 0.25]])
    Pos3 = np.zeros((4,2))
    for j in range (4):
        a = int(Con[0][j])
        x = Pos[a,0]
        y = Pos[a,1]
        Pos3[j][0]= x
        Pos3[j][1]= y
        
    Jaux = Ben*Matrix([[Pos3[0][0],Pos3[0][1]],[Pos3[1][0],Pos3[1][1]],[Pos3[2][0],\
                    Pos3[2][1]],[Pos3[3][0],Pos3[3][1]]])
    J = np.zeros((2,2))
    for i in range (2):
        J[0][i]=Jaux[0,i]
        J[1][i]=Jaux[1,i]
    Jdet = np.linalg.det(J)
    a = np.zeros_like(Efe)
    for j in range (ele2):
        a[j][:] = Ht*Efe[j][:]*Jdet*W1*W2
    return a
    

#%% Función que permite hacer la cuadratura gaussiana.

def Gausscuad():
    p = 0.8611363116
    q = 0.3399810436
    W1 = 0.6521451549 
    W2 = 0.3478548451
    a,a5 = Intgauss(-p,-p,W2,W2)
    b,b5 = Intgauss(-p,p,W2,W2)
    c,c5 = Intgauss(p,-p,W2,W2)
    d,d5 = Intgauss(p,p,W2,W2)
    a2,a6 = Intgauss(-q,-q,W1,W1)
    b2,b6 = Intgauss(-q,q,W1,W1)
    c2,c6 = Intgauss(q,-q,W1,W1)
    d2,d6 = Intgauss(q,q,W1,W1)
    a3,a7 = Intgauss(-p,-q,W2,W1)
    b3,b7 = Intgauss(-p,q,W2,W1)
    c3,c7 = Intgauss(p,-q,W2,W1)
    d3,d7 = Intgauss(p,q,W2,W1)
    a4,a8 = Intgauss(-q,-p,W1,W2)
    b4,b8 = Intgauss(-q,p,W1,W2)
    c4,c8 = Intgauss(q,-p,W1,W2)
    d4,d8 = Intgauss(q,p,W1,W2)

    KeG = a+b+c+d+a2+b2+c2+d2+a3+b3+c3+d3+a4+b4+c4+d4
    MeG = a5+b5+c5+d5+a6+b6+c6+d6+a7+b7+c7+d7+a8+b8+c8+d8
               
    p = 0.8611363116
    q = 0.3399810436
    W1 = 0.6521451549 
    W2 = 0.3478548451
    a = gaussbe(-p,-p,W2,W2)
    b = gaussbe(-p,p,W2,W2)
    c = gaussbe(p,-p,W2,W2)
    d = gaussbe(p,p,W2,W2)
    a2 = gaussbe(-q,-q,W1,W1)
    b2 = gaussbe(-q,q,W1,W1)
    c2 = gaussbe(q,-q,W1,W1)
    d2 = gaussbe(q,q,W1,W1)
    a3 = gaussbe(-p,-q,W2,W1)
    b3 = gaussbe(-p,q,W2,W1)
    c3 = gaussbe(p,-q,W2,W1)
    d3 = gaussbe(p,q,W2,W1)
    a4 = gaussbe(-q,-p,W1,W2)
    b4 = gaussbe(-q,p,W1,W2)
    c4 = gaussbe(q,-p,W1,W2)
    d4 = gaussbe(q,p,W1,W2)
    
    beG = a+b+c+d+a2+b2+c2+d2+a3+b3+c3+d3+a4+b4+c4+d4
    return KeG, MeG, beG


#%% Función que permite hacer la cuadratura gaussiana para elementos triangulares.
#    Ensambla la matriz de rigidez global

def Gausstri():
    e = 0
    n = 0
    W = 2
    
    H = Matrix([[1-e-n, e, n]])
    Ben = Matrix([[-1.,  1.,  0.], [-1.,  0.,  1.]])
    K = np.zeros((nodos,nodos))
    M = np.zeros_like(K)
    be = np.zeros_like(Efe2)
    Pos3 = np.zeros((3,2))
    for k in range (ele):    
        for j in range (3):
                a = int(Con[k][j])
                x = Pos[a,0]
                y = Pos[a,1]
                Pos3[j][0]= x
                Pos3[j][1]= y
        Jaux = Ben*Matrix([[Pos3[0][0],Pos3[0][1]],[Pos3[1][0],Pos3[1][1]]\
                                   ,[Pos3[2][0],Pos3[2][1]]])
        
        J = np.zeros((2,2))
        for i in range (2):
            J[0][i]=Jaux[0,i]
            J[1][i]=Jaux[1,i]
        Jdet = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        B = Jinv*Ben
        Bt = B.transpose()
        Ht = H.transpose()
        Ke = Bt*B*Jdet*W**2
        Me = Ht*H*Jdet*W**2
        for j in range (3):
            for i in range (3):
                a = Ke[j,i]
                y = int(Con[k][i])
                x = int(Con[k][j])
                K[x][y] = a + K[x][y]
                a = Me[j,i]*Ka
                M[x][y] = a + M[x][y]
    
        Ht = [1-e-n, e, n]
        j = 0
        for j in range (3):           
            be[k][j] = Ht[j]*Efe2[k][j]
        
    global Gaux
    Gaux = 1
    return K, M, be


#%% Esta parte se utiliza si la malla es con elementos cuadrados.

#Con2 = cells['quad']
#ele2 = int(np.size(Con2)/4)    

#Efe = np.zeros((ele,4))
#for i in range (ele):
#    for j in range (4):
#        a = int(Con[i][j])
#        x = Pos[a,0]
#        y = Pos[a,1]
#        Efe[i][j] = f(x,y)


#Ke, Me, be = Gausscuad()
#K = np.zeros((nodos,nodos))
#aux = -1
#for k in range (ele):
#    for i in range (4):
#        for j in range (4):
#            aux = aux+1
#            a = Ke[aux,j]
#            y = int(Con[k][j])
#            x = int(Con[k][i])
#            K[x][y] = a + K[x][y]
#            
#
#M = np.zeros_like(K)
#for k in range (ele):
#    for i in range (4):
#        for j in range (4):
#            a = Me[i,j]*Ka
#            y = int(Con[k][j])
#            x = int(Con[k][i])
#            M[x][y] = a + M[x][y]

            
#b = np.zeros((nodos,1))
#for k in range (ele):
#    for i in range (4):
#        a = be[k][i]
#        y = int(Con[k][i])
#        b[y][0] = a 


#%% Esta parte se utiliza si la malla es de elementos triangulares.

Con = cells['triangle']
ele = int(np.size(Con)/3)

Efe2 = np.zeros((ele,3))
for i in range (ele):
    for j in range (3):
        a = int(Con[i][j])
        x = Pos[a][0]
        y = Pos[a][1]
        Efe2[i][j] = f(x,y)
        
K, M, be = Gausstri()

c = np.zeros((nodos,1))
b = np.zeros_like(c)
for k in range (ele):
    for i in range (3):
        a = be[k][i]
        y = int(Con[k][i])
        c[y][0] = a


#%% En este punto se elige qué nodos se ven afectados por el término fuente.

S = K+M

nodoini = 130
nodofin = 190

for a in range (nodoini,nodofin,1):
    d = np.zeros(3)
    d[:] = Con[a][:]
    S[int(d[0])][:] = 0
    S[int(d[0])][int(d[0])] = 1
    b[int(d[0])][0] = c[int(d[0])][0]
    S[int(d[1])][:] = 0
    S[int(d[1])][int(d[1])] = 1
    b[int(d[1])][0] = c[int(d[1])][0]
    S[int(d[2])][:] = 0
    S[int(d[2])][int(d[2])] = 1
    b[int(d[2])][0] = c[int(d[2])][0]
    

#%% Se imponen las condiciones de frontera (dadas desde la información del archivo .msh)

for i in range(np.size(condi)):
    a = condi[i]
    S[a][:] = 0
    S[a][a] = 1
    b[a][0] = 0
    
    
#%% Se resuelve el sistema

u = solve(S,b)


#%% Visualizacion

X = np.zeros(nodos)
Y = np.zeros_like(X)
for i in range(nodos):
    X[i] = Pos[i][0]
    Y[i] = Pos[i][1]
    
plt.scatter(X, Y,linewidth=0.1, c=u[:, 0], cmap="bwr")
plt.show()

fig2 = plt.figure() 
ax = Axes3D(fig2)   
ax.plot_trisurf(X,Y,u[:,0],triangles=Con, cmap='bwr') 

plt.show()



