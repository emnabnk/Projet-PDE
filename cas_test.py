import numpy as np
import matplotlib.pyplot as plt


# Paramètres
d = 0.05# Coefficient de diffusion
c = 0.01# Coefficient de convection

#bornes intervalle
a = -2
b = 2

#ajustement=0.95  
# Pas définis par l'utilisateur
hx = 0.01# Pas en espace
#ht = 0.5*(hx**2/d)*ajustement  # Pas en temps

# Condition de stabilité
gamma=0.4
nu=0.5
delta_t1 = gamma * (hx**2) / d 
delta_t2= 100#nu * hx /c
ht= min(np.abs(delta_t1),np.abs(delta_t2))
print(ht)
print("Ceci est un coucou")
#print(delta_t1)
#print(delta_t2)

# Définition des bornes
Nx = b-a # Longueur de l'intervalle spatial
Nt = 1# Durée de l'intervalle temporel

# Calcul du nombre de points à partir des pas
nbrx = int(Nx / hx) + 1  # Nombre de points en espace
nbrt = int(Nt / ht) + 1  # Nombre de points en temps

# Discrétisation des axes
x = np.linspace(a, b, nbrx)  # Vecteur spatial
t = np.linspace(0, Nt, nbrt)  # Vecteur temporel

# Initialisation des solutions
solutionsApprochées = np.zeros((nbrx, nbrt))
solutionExacte = np.zeros((nbrx, nbrt))

# Fonction solution exacte
def u(x, t):
    return np.sin(x) * np.exp(2*t)

def f(x,t) :
    return 2*np.sin(x)*np.exp(2*t) - d*np.cos(x)*np.exp(2*t) - c*np.sin(x)*np.exp(2*t)

# Conditions initiales
for i in range(nbrx):
    solutionsApprochées[i, 0] = u(x[i], 0)  # Condition initiale en t=0

for j in range(nbrt):
    solutionsApprochées[0, j] = u(a, t[j])  # Condition aux bords en x=0
    solutionsApprochées[nbrx-1, j] = u(b, t[j])   # Condition aux bords en x=1

# Résolution avec différences finies
for n in range(nbrt - 1):  # Parcours temporel
    for i in range(1, nbrx - 1):  # Parcours spatial
        if(c<=0):
            solutionsApprochées[i, n + 1] =solutionsApprochées[i, n] + ht * (d * (solutionsApprochées[i + 1, n] - 2 * solutionsApprochées[i, n] + solutionsApprochées[i - 1, n]) / hx**2 - c * (solutionsApprochées[i + 1, n] - solutionsApprochées[i, n]) / hx + f(x[i],t[n]))
        else :
            solutionsApprochées[i, n + 1] =solutionsApprochées[i, n] + ht * (d * (solutionsApprochées[i + 1, n] - 2 * solutionsApprochées[i, n] + solutionsApprochées[i - 1, n]) / hx**2 - c * (solutionsApprochées[i , n] - solutionsApprochées[i-1, n]) / hx + f(x[i],t[n]))

# Calcul de la solution exacte
for n in range(nbrt):
    for i in range(nbrx):
        solutionExacte[i, n] = u(x[i], t[n])

erreurs2=[]
# Calcul de l'écart

for i in range(nbrx): 
    erreur = np.abs(solutionsApprochées[i,:] - solutionExacte[i,:])
    erreurs2.append

# erreur = np.abs(solutionsApprochées[i,:] - solutionExacte[i,:])
erreur2 = np.linalg.norm(solutionsApprochées - solutionExacte)/np.linalg.norm( solutionExacte)
# # print(f"erreur :{erreur}")
print(f"erreur2:{erreur2}")

# Calcul de l'erreur en fonction du temps
erreurs = []
erreurs1 = []
for n in range(nbrt):
    erreur_t = np.linalg.norm(solutionsApprochées[:,n] - solutionExacte[:,n]) / np.linalg.norm(solutionExacte[:,n])
    erreurs.append(erreur_t)

erreur3 = np.max(erreurs)   
print(f"erreur3:{erreur3}")  



# Affichage des solutions pour différents instants de temps
def tracer_toutes_les_courbes():
    plt.figure(figsize=(12, 8))

    # Choisir des temps spécifiques pour les tracés (5 instants uniformément répartis)
    indices_temps = np.linspace(0, nbrt - 1, 5, dtype=int)

    for n in indices_temps:
        plt.plot(
            x,
            solutionExacte[:, n],
            label=f"Solution exacte (t={t[n]:.2f})",
            linewidth=2
        )
        plt.plot(
            x,
            solutionsApprochées[:, n],
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
tracer_toutes_les_courbes()




# Tracé de l'erreur en fonction du temps
plt.figure(figsize=(10, 6))
plt.loglog(t, erreurs, marker='o', color='red', label="Erreur en fonction du temps")
plt.xlabel("Temps (t)")
plt.ylabel("Erreur relative")
plt.title("Erreur relative entre solution approchée et exacte en fonction du temps")
plt.grid(True)
plt.legend()
plt.show()

erreur_spatiale = np.abs(solutionsApprochées - solutionExacte)
# Tracé de l'écart (carte de chaleur)
plt.figure(figsize=(10, 6))
plt.imshow(erreur_spatiale, extent=[0, Nt, 0, Nx], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label='Écart absolu')
plt.title("Écart entre la solution approchée et la solution exacte")
plt.xlabel("Temps (t)")
plt.ylabel("Espace (x)")
plt.show()