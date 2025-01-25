import numpy as np
import matplotlib.pyplot as plt

######### ETUDE DE LA FONCTION GAUSSIENNE : u_excate = (1+t)*exp(-x**2) ##########

# Initialisation de tous les paramètres

D = 0.05 # Coefficient de diffusion
C = 0 # Coefficient de convection

#intervalle
a=-2
b=2

delat_x = 0.1  # Pas en espace

# Condition de stabilité
gamma=0.1
nu=0.5
delta_t1 = gamma * (delat_x**2) / D 
if (c!=0):
    delta_t2= nu * delta_x /C
else:
    delta_t2 = 10000000
#on recupere le plus petit delta_t
delat_t=min(np.abs(delta_t1),np.abs(delta_t2))
#print(f"delta_t : {delta_t}")

# intervalles saptial et temporel
intervalle_x = b-a 
intervalle_t = 1

# Calcul du nombre de points avec les pas
nbrx = int(intervalle_x / delat_x) + 1  
nbrt = int(intervalle_t / delat_t) + 1  

# Discrétisation spatiale et temporelle
x = np.linspace(a, b, nbrx)  
t = np.linspace(0, intervalle_t, nbrt)  

# Initialisation solutions
sol_app = np.zeros((nbrx, nbrt))
u_exacte = np.zeros((nbrx, nbrt))

#solution exacte
def u(x, t):
    return (1+t) * np.exp(-x**2)

#calcul de la source f
def f(x,t) :
    return np.exp(-x**2) + D*2*x*(1+t)*np.exp(-x**2) - 2*C*(1+t)*np.exp(-x**2)*(1+2*x**2)

#Condition initiale en t=0
for i in range(nbrx):
    sol_app[i, 0] = u(x[i], 0)  

#Conditions aux limites
for j in range(nbrt):
    sol_app[0, j] = u(a, t[j])  #en x=a
    sol_app[nbrx-1, j] = u(b, t[j])   # en x=b

# Calcul de la solution exacte
for n in range(nbrt):
    for i in range(nbrx):
        u_exacte[i, n] = u(x[i], t[n])

# Résolution avec différences finies
for n in range(nbrt - 1):  
    for i in range(1, nbrx - 1):  
        #disjonction de cas en fonction de la valeur de C
        if(C<=0):
            sol_app[i, n + 1] =sol_app[i, n] + delat_t * (D * (sol_app[i + 1, n] - 2 * sol_app[i, n] + sol_app[i - 1, n]) / delat_x**2 - C * (sol_app[i + 1, n] - sol_app[i, n]) / delat_x + f(x[i],t[n]))
        else :
            sol_app[i, n + 1] =sol_app[i, n] + delat_t * (D * (sol_app[i + 1, n] - 2 * sol_app[i, n] + sol_app[i - 1, n]) / delat_x**2 - C * (sol_app[i , n] - sol_app[i-1, n]) / delat_x + f(x[i],t[n]))



#calcul de l'erreur globale
erreur_globale = np.linalg.norm(sol_app - u_exacte)/np.linalg.norm( u_exacte)
print(f"erreur_globale:{erreur_globale}")

# Calcul de l'erreur en fonction du temps
erreurs = []
for n in range(nbrt):
    erreur_t = np.linalg.norm(sol_app[:,n] - u_exacte[:,n]) / np.linalg.norm(u_exacte[:,n])
    erreurs.append(erreur_t)

erreur3 = np.max(erreurs)   
print(f"erreur3:{erreur3}")  

# Calcul de l'erreur en fonction de x
erreurs2=[]
for i in range(nbrx): 
    erreur = np.abs(sol_app[i,:] - u_exacte[i,:])
    erreurs2.append


# Affichage des solutions pour différents instants de temps
def courbes():
    plt.figure(figsize=(12, 8))

    # Choisir des temps spécifiques pour les tracés (5 instants uniformément répartis)
    indices_temps = np.linspace(0, nbrt - 1, 5, dtype=int)

    for n in indices_temps:
        plt.plot(
            x,
            u_exacte[:, n],
            label=f"Solution exacte (t={t[n]:.2f})",
            linewidth=2
        )
        plt.plot(
            x,
            sol_app[:, n],
            linestyle="--",
            label=f"Solution approchée (t={t[n]:.2f})",
            linewidth=2
        )

    plt.title("Solutions exacte et approchée pour différents instants de temps")
    plt.xlabel("Espace (x)")
    plt.ylabel("Valeur de la solution")
    plt.legend()
    plt.grid(True)
    plt.show()

# Tracé des courbes
courbes()


# Tracé de l'erreur en fonction du temps
plt.figure(figsize=(10, 6))
plt.loglog(t, erreurs, marker='o', color='red', label="Erreur en fonction du temps")
plt.xlabel("Temps (t)")
plt.ylabel("Erreur relative")
plt.title("Erreur relative entre solution approchée et exacte en fonction du temps")
plt.grid(True)
plt.legend()
plt.show()

erreur_spatiale = np.abs(sol_app - u_exacte)
# Tracé de l'écart (carte de chaleur)
plt.figure(figsize=(10, 6))
plt.imshow(erreur_spatiale, extent=[0, intervalle_t, 0, intervalle_x], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label='Écart absolu')
plt.title("Écart entre la solution approchée et la solution exacte")
plt.xlabel("Temps (t)")
plt.ylabel("Espace (x)")
plt.show()
