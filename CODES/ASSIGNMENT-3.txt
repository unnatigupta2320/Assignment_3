import numpy as np
 
 
def dir_vec(A,B):
  return B-A
 
def norm_vec(A,B):
  return np.matmul(omat, dir_vec(A,B))
 
#Generate line points
#def line_gen(A,B):
#  len =10
#  dim = A.shape[0]
#  x_AB = np.zeros((dim,len))
#  lam_1 = np.linspace(0,1,len)
#  for i in range(len):
#    temp1 = A + lam_1[i]*(B-A)
#    x_AB[:,i]= temp1.T
#  return x_AB
 
#Generate line intercepts
def line_icepts(n,c):
  e1 = np.array([1,0]) 
  e2 = np.array([0,1]) 
  A = c*e1/(n@e1)
  B = c*e2/(n@e2)
  return A,B
 
#Generate line points
def line_dir_pt(m,A,k1,k2):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(k1,k2,len)
  for i in range(len):
    temp1 = A + lam_1[i]*m
    x_AB[:,i]= temp1.T
  return x_AB
#Generate line points
 
def line_norm_eq(n,c,k):
  len =10
  dim = n.shape[0]
  m = omat@n
  m = m/np.linalg.norm(m)
#  x_AB = np.zeros((dim,2*len))
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(k[0],k[1],len)
#  print(lam_1)
#  lam_2 = np.linspace(0,k2,len)
  if c==0:
    for i in range(len):
      temp1 = lam_1[i]*m
      x_AB[:,i]= temp1.T
  else:
    A,B = line_icepts(n,c)
    for i in range(len):
      temp1 = A + lam_1[i]*m
      x_AB[:,i]= temp1.T
#    temp2 = B + lam_2[i]*m
#    x_AB[:,i+len]= temp2.T
  return x_AB
 
#def line_dir_pt(m,A, dim):
#  len = 10
#  dim = A.shape[0]
#  x_AB = np.zeros((dim,len))
#  lam_1 = np.linspace(0,10,len)
#  for i in range(len):
#    temp1 = A + lam_1[i]*m
#    x_AB[:,i]= temp1.T
#  return x_AB
 
 
#Generate line points
def line_gen(A,B):
  len =10
  x_AB = np.zeros((2,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
 
#Foot of the Altitude
def alt_foot(A,B,C):
  m = B-C
  n = np.matmul(omat,m) 
  N=np.vstack((m,n))
  p = np.zeros(2)
  p[0] = m@A 
  p[1] = n@B
  #Intersection
  P=np.linalg.inv(N.T)@p
  return P
 
#Intersection of two lines
def line_intersect(n1,c1,n2,c2):
  N=np.vstack((n1,n2))
  p = np.array([c1,c2]) 
  #Intersection
  P=np.linalg.inv(N)@p
#  P=np.linalg.inv(N.T)@p
  return P
 
#Radius and centre of the circumcircle
#of triangle ABC
def ccircle(A,B,C):
  p = np.zeros(2)
  n1 = dir_vec(B,A)
  p[0] = 0.5*(np.linalg.norm(A)**2-np.linalg.norm(B)**2)
  n2 = dir_vec(C,B)
  p[1] = 0.5*(np.linalg.norm(B)**2-np.linalg.norm(C)**2)
  #Intersection
  N=np.vstack((n1,n2))
  O=np.linalg.inv(N)@p
  r = np.linalg.norm(A -O)
  return O,r
 
#Radius and centre of the incircle
#of triangle ABC
def icentre(A,B,C,k1,k2):
  p = np.zeros(2)
  t = norm_vec(B,C)
  n1 = t/np.linalg.norm(t)
  t = norm_vec(C,A)
  n2 = t/np.linalg.norm(t)
  t = norm_vec(A,B)
  n3 = t/np.linalg.norm(t)
  p[0] = n1@B- k1*n2@C
  p[1] = n2@C- k2*n3@A
  N=np.vstack((n1-k1*n2,n2-k2*n3))
  I=np.matmul(np.linalg.inv(N),p)
  r = n1@(I-B)
  #Intersection
  return I,r
 
def mult_line(A_I,b_z,k,m):
 for i in range(m):
  if i == 0:
    x = line_norm_eq(A_I[i,:],b_z[i],k[i,:])
  elif i == 1:
    y = line_norm_eq(A_I[i,:],b_z[i],k[i,:])
    z = np.vstack((x[None], y[None]))
  else:
    x = line_norm_eq(A_I[i,:],b_z[i],k[i,:])
    z = np.vstack((z,x[None]))
 return z
 
dvec = np.array([-1,1]) 
#Orthogonal matrix
omat = np.array([[0,1],[-1,0]])

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
 
import sys                                          #for path to external scripts
sys.path.insert(0, '/storage/emulated/0/tlc/school/ncert/linman/codes/CoordGeo') 
 
 
#centre and radius of circles
A=np.array([0,0])
B=np.array([8,0])
r1=4
r2=3

#Input parameters
A = np.array(([0,0]))
B = np.array(([8,0]))
P = np.array(([2,3.4]))
Q= np.array(([2,-3.4]))
R = np.array(([6.875,2.78]))
S = np.array(([6.875,-2.78]))
f1 = -16
f2 = 9
m = np.array(([1,-0.576]))
n = np.array(([1,0.576]))
o=np.array(([1,0.404])) 
p=np.array(([1,-0.404]))
k1 = 8
k2 = -8
l1 = 9
l2 = -9

##Generating all lines
x_AR = line_dir_pt(o,A,l1,l2)
x_AS = line_dir_pt(p,A,l1,l2)
x_BP = line_dir_pt(m,B,k1,k2)
x_BQ = line_dir_pt(n,B,k1,k2)
x_AP = line_gen(A,P)
x_AQ = line_gen(A,Q)
x_BR = line_gen(B,R)
x_BS = line_gen(B,S)
x_AB = line_gen(A,B)
 
 
#Plotting all lines
plt.plot(x_AR[0,:],x_AR[1,:],label='$AR$') 
plt.plot(x_AS[0,:],x_AS[1,:],label='$AS$')
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BP[0,:],x_BP[1,:],label='$BP$')
plt.plot(x_BQ[0,:],x_BQ[1,:],label='$BQ$')
plt.plot(x_AP[0,:],x_AP[1,:],label='$AP$')
plt.plot(x_AQ[0,:],x_AQ[1,:],label='$AQ$')
plt.plot(x_BR[0,:],x_BR[1,:],label='$BR$')
plt.plot(x_BS[0,:],x_BS[1,:],label='$BS$')




plt.plot(A[0], A[1], 'o')
plt.text(0,0,'A',weight = "bold")
plt.plot(B[0], B[1], 'o')
plt.text(8,0,'B',weight = "bold")
plt.plot(P[0], P[1], 'o')
plt.text(2,3.4,'P', weight = "bold")
plt.plot(Q[0], Q[1], 'o')
plt.text(2,-3.4,'Q', weight = "bold")
plt.plot(R[0], R[1], 'o')
plt.text(6.875,2.78,'R', weight = "bold")
plt.plot(S[0], S[1], 'o')
plt.text(6.875,-2.78,'S', weight = "bold")


#Plotting the circles
theta = np.linspace(0,2*np.pi,50)
u=np.array([np.cos(theta),np.sin(theta)])

A1=A.reshape(2,1)
B1=B.reshape(2,1)
C1=A1+r1*u
C2=B1+r2*u

plt.plot(C1[0,:],C1[1,:],label='$circle with center A$')
plt.plot(C2[0,:],C2[1,:],label='$circle with center B$')

 
plt.xlabel('$x axis$')
plt.ylabel('$y axis$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
# plt.plot(range(0, 12))
plt.figure(figsize=(5,5))
plt.show()
 
