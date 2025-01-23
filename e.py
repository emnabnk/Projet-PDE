import numpy as np
import matplotlib.pyplot as plt

# Paramètres
d = 0.5  # Coefficient de diffusion
c = 0.1  # Coefficient de convection
ajustement = 0.95
hx = 0.01  # Pas en espace
gamma = 0.5
nu = 0.5
delta_t1 = gamma * (hx**2) / d
delta_t2 = nu * hx / c
ht = min(np.abs(delta_t1), np.abs(delta_t2)) * ajustement  # Pas en temps
print("Pas temporel (ht) :", ht)

# Définition des bornes
Nx = 1  # Longueur de l'intervalle spatial
Nt = 3  # Durée de l'intervalle temporel
nbrx = int(Nx / hx) + 1  # Nombre de points en espace
nbrt = int(Nt / ht) + 1  # Nombre de points en temps

# Discrétisation des axes
x = np.linspace(0, Nx, nbrx)  # Vecteur spatial
t = np.linspace(0, Nt, nbrt)  # Vecteur temporel

# Initialisation des solutions
solutionsApprochées = np.zeros((nbrx, nbrt))
solutionExacte = np.zeros((nbrx, nbrt))

# Fonction solution exacte
def u(x, t):
    return np.sin(np.pi * x) * (1 + t)

# Terme source
def f(x, t):
    return np.pi**2 * (t + 1) * d * np.sin(np.pi * x) + np.pi * (t + 1) * c * np.cos(np.pi * x) + np.sin(np.pi * x)

# Conditions initiales et aux bords
for i in range(nbrx):
    solutionsApprochées[i, 0] = u(x[i], 0)  # Condition initiale en t=0
for j in range(nbrt):
    solutionsApprochées[0, j] = u(0, t[j])  # Condition au bord x=0
    solutionsApprochées[-1, j] = u(1, t[j])  # Condition au bord x=1

# Résolution avec différences finies
for n in range(nbrt - 1):  # Parcours temporel
    for i in range(1, nbrx - 1):  # Parcours spatial
        if c <= 0:
            solutionsApprochées[i, n + 1] = solutionsApprochées[i, n] + ht * (
                d * (solutionsApprochées[i + 1, n] - 2 * solutionsApprochées[i, n] + solutionsApprochées[i - 1, n]) / hx**2
                - c * (solutionsApprochées[i + 1, n] - solutionsApprochées[i, n]) / hx
                + f(x[i], t[n])
            )
        else:
            solutionsApprochées[i, n + 1] = solutionsApprochées[i, n] + ht * (
                d * (solutionsApprochées[i + 1, n] - 2 * solutionsApprochées[i, n] + solutionsApprochées[i - 1, n]) / hx**2
                - c * (solutionsApprochées[i, n] - solutionsApprochées[i - 1, n]) / hx
                + f(x[i], t[n])
            )

# Calcul de la solution exacte
for n in range(nbrt):
    for i in range(nbrx):
        solutionExacte[i, n] = u(x[i], t[n])

# Calcul des erreurs
erreur_spatiale = np.abs(solutionsApprochées - solutionExacte)
erreur_globale = np.linalg.norm(solutionsApprochées - solutionExacte) / np.linalg.norm(solutionExacte)
print(f"Erreur globale : {erreur_globale}")

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

# Tracé de l'écart (carte de chaleur)
plt.figure(figsize=(10, 6))
plt.imshow(erreur_spatiale, extent=[0, Nt, 0, Nx], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label='Écart absolu')
plt.title("Écart entre la solution approchée et la solution exacte")
plt.xlabel("Temps (t)")
plt.ylabel("Espace (x)")
plt.show()


erreurs_hx = []
hx_values = [0.1, 0.05, 0.025, 0.01]

for hx in hx_values:
    # Recalculer toutes les étapes pour chaque hx (adapter les calculs initiaux)
    # Vous pouvez encapsuler la résolution dans une fonction pour plus de lisibilité
    pass

    erreurs_hx.append(np.linalg.norm(solutionsApprochées - solutionExacte) / np.linalg.norm(solutionExacte))

plt.figure(figsize=(10, 6))
plt.loglog(hx_values, erreurs_hx, marker='o', label="Erreur en fonction de hx")
plt.title("Étude de convergence spatiale")
plt.xlabel("Pas spatial (hx)")
plt.ylabel("Erreur relative")
plt.grid(True, which="both", linestyle="--")
plt.legend()
plt.show()
