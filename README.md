
# Lehmer Pairs Signature

**Local Extraction Method for the Derivatives of the Hardy Z-Function and Analytical Signature of Lehmer Pairs**

## 📄 Abstract

This repository contains the numerical methods, data, and source code accompanying the research note on the stable extraction of the derivatives of the Hardy $Z$-function $Z(t)$. The method relies on double-window integrals combined with error-absorbing ("sponge") linear systems. Applied to Lehmer pairs (anomalously close consecutive zeros of the Riemann zeta function), this approach achieves a precision ranging from $10^{-14}$ to $10^{-20}$ without requiring a fine sweep of the real axis.

Furthermore, we introduce and validate a local curvature formula $S(\gamma_n)$. A comparative study on 500 normal pairs and several Lehmer pairs at high altitudes demonstrates that the relation $\Delta^2(S(t_1) + S(t_2)) \approx 2$ constitutes a robust and unique analytical signature of Lehmer pairs. We provide a rigorous mathematical explanation of this phenomenon through a local/global decomposition of the curvature, including an explicit logarithmic bound on the residual term.

## 📂 Repository Structure

```text
Lehmer-pairs-signature/
├── paper/                          # Research note (LaTeX + figures)
│   ├── note_lehmer.tex             # Main LaTeX source file
│   └── signature_combinee_lehmer_vs_normal.png  # Main combined plot
── data/                           # Raw data
│   ── groupe_temoin_500_paires.csv  # Control group data (500 pairs)
├── scripts/                        # Reproducible Python scripts
│   ├── compute_signature.py        # Script to compute S(gamma_n)
│   └── plot_signature.py           # Script to generate the combined plot
── requirements.txt                # Python dependencies
├── LICENSE                         # MIT License (Code)
── README.md                       # This file
