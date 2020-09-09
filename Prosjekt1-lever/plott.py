import numpy as np
import matplotlib.pyplot as plt


#Leser filen
f= open("v_og_x.txt","r")
list= []
for line in f:
    line1 = float(line[2:])
    list.append(line1)#list[line] = line1



#Deler opp listene
lengde = int((len(list))/3)

x_liste = list[:lengde]
x_liste.insert(0,0)
x_liste.append(1)

v_liste = list[lengde:(2*lengde)]
v_liste.insert(0,0)
v_liste.append(0)

u_liste = list[2*lengde:]
u_liste.insert(0,0)
u_liste.append(0)

#Gjør om listene til arrays
array_list_x = np.array(x_liste)
array_list_v = np.array(v_liste)
array_list_u = np.array(u_liste)


#Plotter den numeriske og matematiske løsningen
print(lengde)
plt.title("grid points n ="+str(lengde))
plt.plot(array_list_x, array_list_v, label="numerical solution")
plt.plot(array_list_x, array_list_u, label="mathematical solution")
plt.xlabel("x")
plt.ylabel("solution (u(x) and v(x))")
plt.legend()
plt.show()
