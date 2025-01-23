# Tracé en log-log des données fournies dans l'image
# POUR C=0.1 et delta_t=0.00002  fixé 

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# Données
hx = np.array([0.05, 0.01, 0.005 ])
erreur = np.array([0.003743, 0.0007598, 0.0004495])
#erreurs_modif=np.array([0.1650, 0.0825, 0.04125, 0.020625,0.0103125 ])

slope, _, _, _, _ = linregress(np.log(hx), np.log(erreur))

print(f"Pente  : {slope}")

# Tracé en log-log
plt.figure(figsize=(10, 6))
plt.loglog(hx, erreur, marker='o', color='blue', label="Erreur en fonction de hx (log-log)")
plt.xlabel("Pas spatial hx (log)")
plt.ylabel("Erreur (log)")
plt.title("Erreur relative en fonction de hx (log-log)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.show()

