from matplotlib import *
from numpy import * 
from sympy import *
from pylab import *
from os import listdir 

theta=[0.0 for i in range(1000)]
A=[0.0 for i in range(1000)]
phi=[0.0 for i in range(1000)]
nombre=["0"for i in range(1000)]

j=0
i=0
n=0
p=0
q=0
l=0
lin=""
#Transpuesta de tvector
def Transponer(D ,m ,n ):
   transpuesta = [[0.0 for j in range(m)] for i in range(n)]
   for i in range(n):
      for j in range(m):   
          transpuesta[i][j] =D[j][i]
   return transpuesta  

#Multiplicacion matriz x matrix
def MultiplicarM(A, filA, colA, V, filV, colV):
  a1 = [[0.0 for m in range(filA)] for n in range(colV)] 
  for i in range(filA):
      for j in range(colV): 
        for k in range(colA):
           a1[i][j]= a1[i][j]+A[i][k]*V[k][j] 
  return a1         

#Multiplicacion matriz X vector        
def MultiplicarV(A, filA, colA, V, filV): 
  b1 = [0.0 for m in range(filA)]  
  for i in range(filA):
     for j in range(colA): 
         b1[i]= b1[i]+A[i][j]*V[j]
  return b1        

#Metodo de Cholesky
def Cholesky(a, filA):
  N = filA
  L1 = [[0.0] * N for i in range(N)]
  for i in range(N):
    acumula = 0.0
    for k in range(i+1):
        acumula = sum(L1[i][j] * L1[k][j] for j in range(k))
        if (i == k): 
           L1[i][k] = sqrt(a[i][i] - acumula)
        else:
           L1[i][k] = ((1.0 / L1[k][k]) * (a[i][k] - acumula)) 
  return L1
                
 # encontrar Y de la ecuacion LY=b
def TrianguloInf(b, filb, L):
  Y1 = [0.0 for j in range(filb)] 
  Y1[0] = b[0]/L[0][0]
  for i in range(1,filb):
    sumar = 0.0
    for j in range(i): 
        sumar = sumar + L[i][j]*Y1[j]
    Y1[i]=(b[i]-sumar)/L[i][i]
  return Y1
       
# encontrar X de la ecuacion LT X =Y
def TrianguloSup(Y, fily, LT):
  X1 = [0.0 for j in range(fily)] 
  X1[fily-1] = Y[fily-1]/LT[fily-1][fily-1]
  for i in range(fily-2, -1, -1):
      sumar = 0.0
      for j in range(fily-1, i, -1): 
          sumar = sumar + LT[i][j]*X1[j]
      X1[i]=(Y[i]-sumar)/LT[i][i]      
  return X1       
##########   PROGRAMA PRINCIPAL  #######################       
        
#leer la lista de archivos
camino= raw_input('\tRuta de la carpeta Brahe-3141-f: ')
path = str(camino)+ "/Brahe-3141-f"
for archivo in listdir(path): 
   a=archivo.split("_")
   nombre[i]=archivo 
   theta[i] = double(a[2])
   phi[i] = double(a[4].replace(".dat", ""))
   i=i+1 
    
#abrir tabla escritura
archi1=open("ReasultadosPunto1.txt","w") 
archi1.write("theta" +"    " +"phi"+"          "+ "y_ini"+"             " +"v_ini"+"              "+ "g"+ "\n\n")

for I in range(1000):
   archi=open(path+"/"+nombre[I].strip(),'r')
   tvector = [[0.0 for j in range(3)] for k in range(38)]
   yvector = [0.0 for k in range(38)] 
   T = [[0.0 for j in range(38)] for k in range(3)] 
   a = [[0.0 for m in range(3)] for n in range(3)]
   TL = [[0.0 for m in range(3)] for n in range(3)] 
   L = [[0.0 for m in range(3)] for n in range(3)] 
   b = [0.0 for m in range(3)]
   Y = [0.0 for j in range(3)] 
   X = [0.0 for j in range(3)]
   k=0  
   for linea in archi.readlines():
       c=linea.split(" ")
       t=double(c[0].replace("\n",""))
       y= double(c[1])
       yvector[k]=y
       tvector[k][0]=1
       tvector[k][1]=t
       tvector[k][2]=t*t/2
       k=k+1
   
   archi.close() 
   T = Transponer(tvector,k,3) 
   a = MultiplicarM(T,3,k, tvector, k,3) #Multiplicacion de transpuesta * tvector
   b = MultiplicarV(T,3,k, yvector, k) #Multiplicacion de transpuesta * yvector 
   L = Cholesky(a, 3)  #Cholesky
   TL = Transponer(L,3,3) #transpuesta L    
   Y = TrianguloInf(b,3,L) # encontrar Y de la ecuacion LY=b
   X = TrianguloSup(Y,3,TL) # encontrar X de la ecuacion LT * X= Y      
   
   lin=lin+str(theta[I]) +"    " +str(phi[I])+"    "+ str(X[0])+"    " +str( X[1])+"    "+ str( X[2])+ "\n"


   if(X[0]<2.0):  
    q = q + 1
    A[I]=X[2]
    
   #print theta[I],"  " ,phi[I]," ",  X[0],"  " ,  X[1],"  ", X[2]    

archi1.write(lin)
#archi1.close()
    
############################################Segundo punto###########################################################

G=[0.0 for i in range(q)]
theta1=[0.0 for i in range(q)]
for i in range (1000):     #se descartan los datos con y  mayores que 2
 if A[i]!=0.0:
    G[l]=A[i]
    theta1[l]=theta[i]
    l=l+1


fig1=scatter(theta1,G)
xlabel('theta')
ylabel('gravedad')
title(u'g vs theta')
savefig('GraficaPunto2.png')
close()
#########################################Tercer punto###########################################################


F=[0.0 for i in range(q)]
suma=(sum(G[i] for i in range(q)))/q
for i in range (q):
   F[i]=1-G[i]/suma

#imprimir resultados             
lin=""
for i in range (q):
  lin=lin+str(F[i])+"      "+str(theta1[i])+"\n"

archi2=open("ReasultadosPunto3.txt","w") 
archi2.write(" variaciones f" +"        " +"thata"+ "\n\n") 
archi2.write(lin)


#graficar resultados
fig2=scatter(theta1,F)
xlabel('theta')
ylabel('variaciones de la gravedad')
title(u'F vs theta')
savefig('GraficaPunto3.png')
close()        
#############################################Cuarto punto########################################################

###primera ecuacion #####
matriz1 = [[0.0 for i in range(2)] for j in range(q)]  
for k in range(q):  

  matriz1[k][0]=1
  matriz1[k][1]=cos(2*theta1[k])
 
T = Transponer(matriz1,q,2) 
a = MultiplicarM(T,2,q, matriz1, q,2) #Multiplicacion de transpuesta * tvector
b = MultiplicarV(T,2,q, F, q) #Multiplicacion de transpuesta * yvector 
L = Cholesky(a, 2)  #Cholesky
TL = Transponer(L,2,2) #transpuesta L    
Y = TrianguloInf(b,2,L) # encontrar Y de la ecuacion LY=b
X = TrianguloSup(Y,2,TL) # encontrar X de la ecuacion LT * X= Y"""

#imprimir resultados
lin=""
lin="a0="+str(X[0])+"  " +"a1="+str(X[1])+"\n"
archi3=open("ReasultadosPunto4.txt","w") 
archi3.write(lin)


###segunda ecuacion####

matriz2 = [[0.0 for i in range(3)] for j in range(q)]  
for k in range(q):  

  matriz2[k][0]=1
  matriz2[k][1]=theta1[k]
  matriz2[k][2]=(theta1[k])**2

T = Transponer(matriz2,q,3) 
a = MultiplicarM(T,3,q, matriz2, q,3) #Multiplicacion de transpuesta * tvector
b = MultiplicarV(T,3,q, F, q) #Multiplicacion de transpuesta * yvector 
L = Cholesky(a, 3)  #Cholesky
TL = Transponer(L,3,3) #transpuesta L    
Y = TrianguloInf(b,3,L) # encontrar Y de la ecuacion LY=b
X = TrianguloSup(Y,3,TL) # encontrar X de la ecuacion LT * X= Y"""

#imprimir resultados
lin=""
lin="b0="+str(X[0])+"  " +"b1="+str(X[1])+"  " +"b2="+str(X[2])+"\n"
archi3.write(lin)


  
  
   
