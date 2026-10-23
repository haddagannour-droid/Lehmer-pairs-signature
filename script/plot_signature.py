"""
plot_signature.py
Génère le graphique combiné : groupe témoin (500 paires) vs paires de Lehmer validées.
Usage : python plot_signature.py
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path  # <--- AJOUT IMPORTANT

# ============================================================
# 1. CHARGEMENT DES DONNÉES (Chemin robuste)
# ============================================================
# On trouve le dossier où se trouve ce script, puis on remonte d'un cran vers 'data'
script_dir = Path(__file__).parent.resolve()
csv_path = script_dir.parent / 'data' / 'groupe_temoin_500_paires.csv'

print(f"Recherche du fichier CSV dans : {csv_path}")

if not csv_path.exists():
    raise FileNotFoundError(f"Le fichier CSV est introuvable à l'adresse : {csv_path}\n"
                            "Vérifiez que le fichier 'groupe_temoin_500_paires.csv' est bien dans le dossier 'data/'.")

df = pd.read_csv(csv_path)
df_normal = df[df['Est_Normale'] == True]

# ============================================================
# 2. PAIRES DE LEHMER VALIDÉES
# ============================================================
lehmer_pairs = [
    {"t1": 7005.062866174921,       "delta2_S": 2.01386, "label": "n°6709-6710"},
    {"t1": 23153.51496722261,       "delta2_S": 2.04399, "label": "n°11410-11411"},
    {"t1": 220538.8529632685,       "delta2_S": 2.00431, "label": "t ~ 2.2e5"},
    {"t1": 947107.8201191042,       "delta2_S": 2.01000, "label": "t ~ 9.5e5"},
    {"t1": 4316832.467101349,       "delta2_S": 2.00027, "label": "t ~ 4.3e6"},
    {"t1": 8436151843.542355,       "delta2_S": 2.00929, "label": "t ~ 8.4e9"},
    {"t1": 7954022502373.43289,     "delta2_S": 2.00000, "label": "t ~ 7.95e12"},
]

# ============================================================
# 3. CONFIGURATION ET TRACÉ
# ============================================================
fig, ax = plt.subplots(figsize=(14, 8))

# Groupe témoin (fond bleu)
ax.scatter(df_normal['t1'], df_normal['Delta2_Somme_S'],
           color='royalblue', alpha=0.4, s=25,
           label=f'Paires normales ({len(df_normal)} paires)')

# Paires de Lehmer (étoiles rouges)
t1_lehmer = [p["t1"] for p in lehmer_pairs]
s_lehmer = [p["delta2_S"] for p in lehmer_pairs]
ax.scatter(t1_lehmer, s_lehmer,
           color='crimson', edgecolor='black', s=180, marker='*',
           zorder=5, label='Paires de Lehmer (Validées)')

# Annotations
for i, txt in enumerate([p["label"] for p in lehmer_pairs]):
    ax.annotate(txt, (t1_lehmer[i], s_lehmer[i]),
                textcoords="offset points", xytext=(10, 5),
                fontsize=9, color='darkred', fontweight='bold')

# Lignes de référence
ax.axhline(y=2.0, color='forestgreen', linestyle='--', linewidth=2.5,
           label='Signature Théorique ($y = 2.0$)', zorder=4)
ax.axhline(y=np.pi**2/3, color='orange', linestyle=':', linewidth=2,
           label=f'Moyenne théorique ($\\pi^2/3 \\approx 3.29$)', zorder=4)

# ============================================================
# 4. MISE EN FORME
# ============================================================
ax.set_xscale('log')
ax.set_xlim(10, 1e14)
ax.set_ylim(-5, 15)
ax.set_xlabel('Altitude $t_1$ (échelle logarithmique)', fontsize=13)
ax.set_ylabel('$\\Delta^2 (S(t_1) + S(t_2))$', fontsize=13)
ax.set_title('Signature Analytique : Convergence exclusive des paires de Lehmer vers 2.0',
             fontsize=15, fontweight='bold', pad=15)
ax.grid(True, which="both", ls="--", alpha=0.5, zorder=0)
ax.legend(loc='upper right', fontsize=11, framealpha=0.9)

# ============================================================
# 5. SAUVEGARDE
# ============================================================
output_path = '../paper/signature_combinee_lehmer_vs_normal.png'
plt.tight_layout()
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"✅ Graphique sauvegardé sous '{output_path}'")
plt.show()