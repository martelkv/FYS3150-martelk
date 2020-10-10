import numpy as np
import matplotlib.pyplot as plt

"""
f= open("textfil.txt","r")
sola,jorda = np.loadtxt(f,usecols=(0,1), unpack=True)


n = int(len(sola)/3)
#print(n)
#print(jorda)
#print(jorda[:n])
#print(jorda[n:2*n])

plt.plot(sola[:n], sola[n:2*n], "o",label="Sun")
plt.plot(jorda[:n], jorda[n:2*n], label="Earth")
plt.legend()
plt.show()
"""
"""
n =2
Numberofstep= 1000
x=np.zeros((Numberofstep,n))
y=np.zeros((Numberofstep,n))
z=np.zeros((Numberofstep,n))

for i in range(n):
    filnavn= "Euler_"+str(i)+".txt"
    x[:,i],y[:,i],z[:,i] = np.loadtxt(filnavn,usecols=(0,1,2), unpack=True)

plt.plot(x[:,0],y[:,0],"o")
plt.plot(x[:,1],y[:,1])


plt.show()
"""
plt.figure(figsize=(6,6))
x,y,z = np.loadtxt("Euler_0.txt", usecols=(0,1,2), unpack =True)
plt.plot(x,y,"o")
x,y,z = np.loadtxt("Euler_1.txt", usecols=(0,1,2), unpack =True)
plt.plot(x,y)

plt.show()
