# -*- coding: utf-8 -*-

"""
Created on Tue May 22 13:51:02 2018

@author: Ana Sofía  Marulanda Duque, Manuel Montoya Zuluaga, Maria Josef Lopera Acosta. 
"""

import numpy as np
from sympy import symbols,Matrix,diff,integrate,sin,cos
from scipy.linalg import solve
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import Tk
from tkinter import Frame, Label, Entry, Button, Scrollbar, END, E, W
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt   #Se importan las librerías que se van a utilizar
from mpl_toolkits.mplot3d import Axes3D
import meshio 
from matplotlib import figure
#%%%Ventana
root= Tk() #Crea la ventana 
root.wm_title("Análisis acústico por elementos finitos")
root.config(background = "GhostWhite")
root.geometry('800x510')    
root.resizable(False, False)  
#%%
#Frame de la izquierda
leftFrame = Frame(root, width=300, height = 700, padx=10)
leftFrame.grid(row=3, column=0, pady=0, padx=15)
#Frame para el título 
Upframe= Frame(root, width= 700, height= 200)
Upframe.grid(row=0, column=0, columnspan=8, padx= 20, pady=12)
titulo=Label (Upframe, text="Propagación del sonido en un teatro", fg="royal blue", font=('Elephant', 20)).grid(row=0, column=0)
subt=Label(Upframe, font= ("Franklin Gothic Book", 11), text="El análisis actústico del problema de la distribución sonora \n en un teatro es solucionado mediante el método de elementos finitos\n que brinda simplicidad, eficacia y precisión.").grid(row=2, column=0, pady=5)

#%%%

#Parámetros del modelo. 
config= Label(leftFrame, font=("Elephant", 13), fg= "royal blue", text="Configuraciones del modelo").grid(row=2, column=0, padx=10, pady=5)


#Scrollbar con las opciones de geometrías. 
Label(leftFrame, text="Geometría:",font=('Franklin Gothic Book', 10)).grid(row=22, column=0, padx=10, pady=2)
RightFrame= Frame(root, width=400, height=500)
RightFrame.grid(row=3, column=1)
grafica= Label(RightFrame, text="Distribución sonora para la geometría", font=('Franklin Gothic Book', 10)).grid(row=0, column=0, sticky=EW)
graf1 = Canvas(RightFrame, width=450, height=320, bg='white').grid(row=1, column=0, padx=10, pady=2,sticky=EW)


#%%
def geo1():
    global Con
    global u 
    global Pos
    points, cells, point_data, cell_data, field_data = meshio.read('teatro1.msh')
    g = cells['line']
    condi = list(set(g.flatten()))
    
    
    
    w = 1
    v = 1
    Ka = w**2/v**2
    
    def f(x,y):
        f = x*y
    #    if (0.48<x<0.52 and 0.3<y<0.32):    
    #        f = 2
    #    else:
    #        f = 0
        return f
         
    
    def Intgauss(e,n,W1,W2):
        H = Matrix([[1/4*(1-e)*(1-n), 1/4*(1+e)*(1-n), 1/4*(1+e)*(1+n), 1/4*(1-e)*(1+n)]])
        Ben = Matrix([[0.25*n - 0.25, -0.25*n + 0.25, 0.25*n + 0.25, -0.25*n - 0.25],\
              [0.25*e - 0.25, -0.25*e - 0.25, 0.25*e + 0.25, -0.25*e + 0.25]])
        Pos3 = np.zeros((4,2))
        Ke = np.zeros((ele*4,4))
        for j in range (4):
            a = int(Con[0][j])
            x = Pos[a,0]
            y = Pos[a,1]
            Pos3[j][0]= x
            Pos3[j][1]= y
        Jaux = Ben*Matrix([[Pos3[0][0],Pos3[0][1]],[Pos3[1][0],Pos3[1][1]],[Pos3[2][0],\
                        Pos3[2][1]],[Pos3[3][0],Pos3[3][1]]])
        
    
        aux = -4
        for i in range(ele):
            aux = aux + 4
            for j in range (4):
                a = int(Con[i][j])
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
            c = a[0,:]
            Ke[aux][:] = c
            Ke[aux+1][:] = a[1,:]
            Ke[aux+2][:] = a[2,:]
            Ke[aux+3][:] = a[3,:]
            
            Me[aux][:] = b[0,:]
            Me[aux+1][:] = b[1,:]
            Me[aux+2][:] = b[2,:]
            Me[aux+3][:] = b[3,:]
            
        return Ke, Me
            
            
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
        for j in range (ele):
            a[j][:] = Ht*Efe[j][:]*Jdet*W1*W2
        return a
        
    
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
    
    
    def Gausstri():
        e = 0
        n = 0
        W = 2
        
        H = Matrix([[1-e-n, e, n]])
        Ben = Matrix([[-1.,  1.,  0.], [-1.,  0.,  1.]])
        K = np.zeros((nodos,nodos))
        M = np.zeros_like(K)
        be = np.zeros_like(Efe2)
    #    Jaux = Ben*Matrix([[Pos[0,0],Pos[0,1]],[Pos[1,0],Pos[1,1]],[Pos[nodosx,0],\
    #                        Pos[nodosx,1]]])
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
            be[k][:] = Ht*Efe2[k][:]*Jdet*W**2
        
        global Gaux
        Gaux = 1
        return K, M, be
    
    
    Con = cells['triangle']
    ele = int(np.size(Con)/3)
    Pos = points
    nodos = int (np.size(Pos)/3)
    
    #h = np.zeros(tam)
    #for i in range(tam):
    #    h[i] = Pos[i][0] 
    #h = list(set(h))
    #hmax = max(h)
    
    
    #Con2 = np.zeros((ele*2,3))
    #aux = -1
    #for j in range (0,ele*2,2):
    #    aux = aux+1
    #    Con2[j][0] = Con[aux][0]
    #    Con2[j][1] = Con[aux][1]
    #    Con2[j][2] = Con[aux][3]
    #    Con2[j+1][0] = Con[aux][1]
    #    Con2[j+1][1] = Con[aux][3]
    #    Con2[j+1][2] = Con[aux][2]
    
    
    #Efe = np.zeros((ele,4))
    #for i in range (ele):
    #    for j in range (4):
    #        a = int(Con[i][j])
    #        x = Pos[0,a]
    #        y = Pos[1,a]
    #        Efe[i][j] = f(x,y)
    
    
    Efe2 = np.zeros((ele,3))
    for i in range (ele):
        for j in range (3):
            a = int(Con[i][j])
            x = Pos[a][0]
            y = Pos[a][1]
            Efe2[i][j] = f(x,y)
            
 
    
    #Ke, Me, be = Gausscuad()
    K, M, be = Gausstri()
     
    
    #
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
    
    b = np.zeros((nodos,1))
    for k in range (ele):
        for i in range (3):
            a = be[k][i]
            y = int(Con[k][i])
            b[y][0] = a
    
    S = K+M
    
    #for a in range (55,70):
    #    S[a][:] = 0
    #    S[a][a] = 1
    #    b[a][0] = f(Pos[a][0],Pos[a][0])
    
    
    
    
    for i in range(np.size(condi)):
        a = condi[i]
        S[a][:] = 0
        S[a][a] = 1
        b[a][0] = 0
        
        
    
    u = solve(S,b)
    
    X = np.zeros(nodos)
    Y = np.zeros_like(X)
    for i in range(nodos):
        X[i] = Pos[i][0]
        Y[i] = Pos[i][1]
    
    
    #umax = max(u)
    #umin = min(u)
    #um = max([abs(umin),abs(umax)])
    #for i in range (nodos):
    #    a = u[i]/um
    #    a = abs(a)
    #    if (a == 0):
    #        a = 0.05
    
    #puntitos=plt.Figure(figsize=(3,2))
    #punti=plt.scatter(X, Y,linewidth=0.1, c=u[:, 0], cmap="bwr")
   
    fig2 = plt.figure(figsize=(4,3)) #Se crea la figura
    ax = Axes3D(fig2)   #Se pone la figura en 3D
    ax.plot_trisurf(X,Y,u[:,0],triangles=Con, cmap='bwr') 
    grafica2= FigureCanvasTkAgg(fig2, master=RightFrame)
    #grafica2.show
    grafica2.get_tk_widget().grid(row=1, column=0)
    
    return u, Con, Pos 
#%% ARCHIVO VTK 
point_data = {'Datos': u}
quad_mesh = {
  'points': Pos,
  'cells': {"triangle": Con}}
meshio.write('Geometria1.vtk', quad_mesh['points'], quad_mesh['cells'], point_data=point_data)
    
#%%
geom1= Button(leftFrame, text="Geometría 1", state= NORMAL, font=('Franklin Gothic Book', 9), command=botongeometria1).grid(row=23, column=0, sticky=EW)
geom1= Button(leftFrame, text="Geometría 2", state= DISABLED, font=('Franklin Gothic Book', 9)).grid(row=24, column=0, sticky=EW)
geom1= Button(leftFrame, text="Geometría 3", state= DISABLED, font=('Franklin Gothic Book', 9)).grid(row=25, column=0, sticky=EW)

ingreseaca=Label(leftFrame, text="Ingresar geometría (gmsh)",font=('Franklin Gothic Book', 10)).grid(row=20, column=0)
geoing=Entry(leftFrame, state= DISABLED)
geoing.insert(END, ".msh")
geoing.grid(row=21, column=0)
def getgeo():
    global geousuario
    geousuario= geoing.get()
    print (geousuario)
    return geousuario
botongeo= Button(leftFrame, text="Listo", command=getgeo, state=NORMAL,font=('Franklin Gothic Book', 9)).grid(row=21, column=0, sticky=E)

mainloop()