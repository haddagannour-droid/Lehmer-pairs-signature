"""
compute_signature.py
Calcule S(γₙ) pour une paire de zéros donnée à partir des dérivées extraites.
Usage : python compute_signature.py
"""

import mpmath as mp

# ============================================================
# PARAMÈTRES
# ============================================================
mp.mp.dps = 100  # Précision à 100 décimales

# Paires de Lehmer validées (c0, t1, t2)
paires_validees = [
    (mp.mpf('7005.081855096608'), mp.mpf('7005.062866174921'), mp.mpf('7005.100564672647')),
    (mp.mpf('23153.54445344124'), mp.mpf('23153.51496722261'), mp.mpf('23153.57422707674')),
    (mp.mpf('220538.8615505244'), mp.mpf('220538.8529632685'), mp.mpf('220538.8702102205')),
    (mp.mpf('947107.8263269367'), mp.mpf('947107.8201191042'), mp.mpf('947107.8325415508')),
    (mp.mpf('4316832.468262171'), mp.mpf('4316832.467101349'), mp.mpf('4316832.469424752')),
    (mp.mpf('8436151843.546886'), mp.mpf('8436151843.542355'), mp.mpf('8436151843.551400')),
    (mp.mpf('7954022502373.432892546996'), mp.mpf('7954022502373.432890153872'), mp.mpf('7954022502373.432894940123')),
]

def fmt(val, n=6):
    return mp.nstr(val, n, strip_zeros=False)

# ============================================================
# FONCTIONS MATHÉMATIQUES
# ============================================================
def C_prime_over_C(t):
    """Dérivée logarithmique C'/C (analytique)"""
    return 2 * t / (t**2 + 0.25) - 0.5 * mp.im(mp.digamma(0.25 + 0.5j * t))

def C_double_prime_over_C(t):
    """Dérivée seconde logarithmique C''/C (analytique)"""
    Cp_over_C = C_prime_over_C(t)
    trigamma_val = mp.psi(1, 0.25 + 0.5j * t)
    d2_ln_C = 2 * (0.25 - t**2) / (t**2 + 0.25)**2 - 0.25 * mp.re(trigamma_val)
    return Cp_over_C**2 + d2_ln_C

def calculer_S(gamma_n, Z1, Z2, Z3):
    """Formule de courbure locale S(γₙ)"""
    Cp_C = C_prime_over_C(gamma_n)
    Cpp_C = C_double_prime_over_C(gamma_n)
    Z2_Z1 = Z2 / Z1
    Z3_Z1 = Z3 / Z1
    return (Cp_C)**2 - 0.5 * Cp_C * Z2_Z1 + 0.25 * (Z2_Z1)**2 \
           - 1.5 * Cpp_C - 0.5 * Z3_Z1 - 1 / (4 * gamma_n**2)

def Z_func(t):
    """Fonction Z de Hardy"""
    th = mp.im(mp.loggamma(0.25 + 0.5j * t)) - (t / 2) * mp.log(mp.pi)
    return mp.re(mp.e**(1j * th) * mp.zeta(0.5 + 1j * t))

# ============================================================
# EXÉCUTION
# ============================================================
print("CALCUL DE S(γ) POUR LES PAIRES DE LEHMER VALIDÉES")
print("=" * 70)

for i, (c0, t1, t2) in enumerate(paires_validees, 1):
    print(f"\n--- Paire {i} (t1 ≈ {fmt(t1, 10)}, Δ = {fmt(t2 - t1, 6)}) ---")
    
    delta = t2 - t1
    
    # Extraction des dérivées via mp.diff
    Z1_t1 = mp.diff(Z_func, t1, 1)
    Z2_t1 = mp.diff(Z_func, t1, 2)
    Z3_t1 = mp.diff(Z_func, t1, 3)
    S_t1 = calculer_S(t1, Z1_t1, Z2_t1, Z3_t1)
    
    Z1_t2 = mp.diff(Z_func, t2, 1)
    Z2_t2 = mp.diff(Z_func, t2, 2)
    Z3_t2 = mp.diff(Z_func, t2, 3)
    S_t2 = calculer_S(t2, Z1_t2, Z2_t2, Z3_t2)
    
    somme_S = S_t1 + S_t2
    test_S = (delta**2) * somme_S
    
    print(f"  S(t1)           = {fmt(S_t1, 6)}")
    print(f"  S(t2)           = {fmt(S_t2, 6)}")
    print(f"  Δ²(S₁+S₂)       = {fmt(test_S, 6)}  (doit être ≈ 2.0)")

print(f"\n[Finished]")