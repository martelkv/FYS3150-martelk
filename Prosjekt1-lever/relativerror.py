import numpy as np
import matplotlib.pyplot as plt


#Leser filen, reger ut h og maksimal error, og legger det til i en liste
h_list=[]
eps_verdi = []
for i in range(6):
    filnavn= "resultat_"+str(i+1)+".txt"
    eps = np.loadtxt(filnavn,usecols=(0), unpack=True)
    eps = np.log10(eps)
    n = len(eps)
    h = 1/(n+1)
    h = np.log10(h)
    h_list.append(h)
    eps_max = np.max(eps)
    eps_verdi.append(eps_max)

#print(eps_verdi)
#print(h_list)

#Plotter h mot maksimal error
plt.plot(h_list, eps_verdi, label= "relative error")
plt.legend()
plt.xlabel("log10(step lenght)")
plt.ylabel("log10(max value of the relative error)")
plt.title("the relative error as function of log10(h) for the function values ui and v")
plt.show()
