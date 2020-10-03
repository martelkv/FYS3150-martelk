import numpy as np
import matplotlib.pyplot as plt
# Dette programmet plotter egenvektoren regnet ut i "Prosjekt2.cpp"

#Leser filen
f= open("plott_egenvektor.txt","r")

#Henter egenvektoren fra filen.
egenvektor = np.loadtxt(f,usecols=(0), unpack=True)

n = len(egenvektor)
print(n)
lengde = np.linspace(0,n-1,n)
print(lengde)
print(egenvektor)


#Plotter egenvektoren mot nummeret på egenvektoren.
plt.plot(lengde,egenvektor)
plt.xlabel("Egenvektornummer")
plt.ylabel("Egenvektor")
plt.title("Formen på stangen for den laveste egenverdien, n = 100")
plt.show()
