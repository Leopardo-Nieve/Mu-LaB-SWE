"""
plot_numerical_results.py
Génère les mêmes contour plots que mms_sympy.py mais pour les résultats numériques.

Lit automatiquement tous les CSV r*_dx*.csv dans ./results/convergence_spatial/
et génère un plot 2x3 (h, u, v, erreur_h, erreur_u, erreur_v) pour chaque niveau.

Usage:
    python plot_numerical_results.py

Sorties dans ./results/convergence_spatial/ :
    - numerical_plot_r{r}_dx{dx}.png  : plot par niveau de raffinement
    - numerical_plot_combined.png      : tous les niveaux côte à côte
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CONV_DIR  = "./results/convergence_spatial"
OUT_DIR   = "./results/convergence_spatial"

COL_X     = "x (m)"
COL_Y     = "y (m)"
COL_H     = "h (m)"
COL_U     = "u (m/s)"
COL_V     = "v (m/s)"
COL_HANAL = "h analytical (m)"
COL_UANAL = "u analytical (m/s)"

# ---------------------------------------------------------------------------
# Détection des CSV
# ---------------------------------------------------------------------------
pattern = re.compile(r"^r(\d+)_dx([\d.]+)\.csv$")
entries = []
for fname in os.listdir(CONV_DIR):
    m = pattern.match(fname)
    if m:
        r_val = int(m.group(1))
        dx_val = float(m.group(2))
        entries.append((r_val, dx_val, os.path.join(CONV_DIR, fname)))

entries.sort(key=lambda e: e[1], reverse=True)  # grossier → fin

if not entries:
    raise RuntimeError(f"Aucun fichier r*_dx*.csv trouvé dans {CONV_DIR}/")

print(f"Fichiers détectés ({len(entries)} niveaux):")
for r, dx, path in entries:
    print(f"  r={r}  dx={dx:.6f}  →  {os.path.basename(path)}")

# ---------------------------------------------------------------------------
# Fonction de plot pour un niveau
# ---------------------------------------------------------------------------
def plot_level(r, dx, csv_path, save_path):
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df.columns = df.columns.str.strip()

    x_vals = df[COL_X].values
    y_vals = df[COL_Y].values
    h_num  = df[COL_H].values
    u_num  = df[COL_U].values
    h_anal = df[COL_HANAL].values
    u_anal = df[COL_UANAL].values

    # v analytique : 0.015*cos(pi*x/Lx)*sin(2*pi*y/Ly)
    Lx_dom = x_vals.max() + (x_vals[1] - x_vals[0]) * 0.5  # approx domainX
    Ly_dom = y_vals.max() + (y_vals[1] - y_vals[0]) * 0.5  # approx domainY
    v_anal = 0.015 * np.cos(np.pi * x_vals / Lx_dom) * np.sin(2 * np.pi * y_vals / Ly_dom)
    v_num  = df[COL_V].values

    err_h = h_num - h_anal
    err_u = u_num - u_anal
    err_v = v_num - v_anal

    # Grille pour contourf — compter valeurs uniques
    Nx = df["x (nodes)"].nunique()
    Ny = df["y (nodes)"].nunique()
    xs = np.linspace(x_vals.min(), x_vals.max(), Nx)
    ys = np.linspace(y_vals.min(), y_vals.max(), Ny)

    def to_grid(vals):
        # CSV écrit x en boucle externe, y en boucle interne
        return vals.reshape(Nx, Ny).T

    H  = to_grid(h_num)
    U  = to_grid(u_num)
    V  = to_grid(v_num)
    Ha = to_grid(h_anal)
    Ua = to_grid(u_anal)
    Va = to_grid(v_anal)
    Eh = to_grid(err_h)
    Eu = to_grid(err_u)
    Ev = to_grid(err_v)

    XX, YY = np.meshgrid(xs, ys)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(
        f"Résultats numériques — r={r}, dx={dx:.5f} m\n"
        r"$h_{num}$, $u_{num}$, $v_{num}$ et erreurs $h_{num}-h_{MMS}$",
        fontsize=12
    )

    def sym_norm(data):
        vmax = max(abs(data.min()), abs(data.max()))
        if vmax == 0:
            vmax = 1e-10
        return TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fields = [
        (H,  r'$h_{num}$ (m)',              'Blues',   False),
        (U,  r'$u_{num}$ (m/s)',             'RdBu_r',  True),
        (V,  r'$v_{num}$ (m/s)',             'RdBu_r',  True),
        (Eh, r'$h_{num} - h_{MMS}$ (m)',     'RdBu_r',  True),
        (Eu, r'$u_{num} - u_{MMS}$ (m/s)',   'RdBu_r',  True),
        (Ev, r'$v_{num} - v_{MMS}$ (m/s)',   'RdBu_r',  True),
    ]

    for ax, (field, title, cmap, diverging) in zip(axes.flat, fields):
        if diverging:
            norm = sym_norm(field)
            cf = ax.contourf(XX, YY, field, 20, cmap=cmap, norm=norm)
        else:
            cf = ax.contourf(XX, YY, field, 20, cmap=cmap)
        plt.colorbar(cf, ax=ax, format='%.4f')
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"  Plot sauvegardé : {save_path}")
    plt.close()

    return {
        'r': r, 'dx': dx,
        'L2_h': np.sqrt(np.mean(err_h**2)),
        'L2_u': np.sqrt(np.mean(err_u**2)),
        'L1_h': np.mean(np.abs(err_h)),
        'L1_u': np.mean(np.abs(err_u)),
        'H': H, 'U': U, 'V': V,
        'Eh': Eh, 'Eu': Eu, 'Ev': Ev,
        'XX': XX, 'YY': YY,
    }

# ---------------------------------------------------------------------------
# Plots individuels
# ---------------------------------------------------------------------------
results = []
for r, dx, path in entries:
    save = os.path.join(OUT_DIR, f"numerical_plot_r{r}_dx{dx:.8f}.png")
    res = plot_level(r, dx, path, save)
    results.append(res)

# ---------------------------------------------------------------------------
# Plot combiné : une ligne par niveau, 3 colonnes (h, u, erreur_h)
# ---------------------------------------------------------------------------
n = len(results)
fig, axes = plt.subplots(n, 3, figsize=(14, 4 * n))
if n == 1:
    axes = [axes]

fig.suptitle("Comparaison numérique / MMS — tous les niveaux de raffinement",
             fontsize=13, fontweight='bold')

for row, res in enumerate(results):
    XX, YY = res['XX'], res['YY']
    dx = res['dx']
    r  = res['r']

    def sym_norm(data):
        vmax = max(abs(data.min()), abs(data.max()))
        if vmax == 0: vmax = 1e-10
        return TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # h numérique
    cf = axes[row][0].contourf(XX, YY, res['H'], 20, cmap='Blues')
    plt.colorbar(cf, ax=axes[row][0], format='%.4f')
    axes[row][0].set_title(f'r={r}, dx={dx:.4f} — $h_{{num}}$ (m)', fontsize=10)
    axes[row][0].set_xlabel('x (m)'); axes[row][0].set_ylabel('y (m)')

    # u numérique
    cf = axes[row][1].contourf(XX, YY, res['U'], 20, cmap='RdBu_r',
                                norm=sym_norm(res['U']))
    plt.colorbar(cf, ax=axes[row][1], format='%.4f')
    axes[row][1].set_title(f'r={r}, dx={dx:.4f} — $u_{{num}}$ (m/s)', fontsize=10)
    axes[row][1].set_xlabel('x (m)'); axes[row][1].set_ylabel('y (m)')

    # erreur h
    cf = axes[row][2].contourf(XX, YY, res['Eh'], 20, cmap='RdBu_r',
                                norm=sym_norm(res['Eh']))
    plt.colorbar(cf, ax=axes[row][2], format='%.2e')
    axes[row][2].set_title(
        f'r={r}, dx={dx:.4f} — $h_{{num}}-h_{{MMS}}$  '
        f'L2={res["L2_h"]:.3e}', fontsize=10)
    axes[row][2].set_xlabel('x (m)'); axes[row][2].set_ylabel('y (m)')

plt.tight_layout()
combined_path = os.path.join(OUT_DIR, "numerical_plot_combined.png")
plt.savefig(combined_path, dpi=150, bbox_inches='tight')
print(f"  Plot combiné sauvegardé : {combined_path}")
plt.close()

# ---------------------------------------------------------------------------
# Résumé erreurs
# ---------------------------------------------------------------------------
print()
print(f"  {'r':>4}  {'dx':>10}  {'L1(h)':>12}  {'L2(h)':>12}  {'L1(u)':>12}  {'L2(u)':>12}")
print(f"  {'-'*58}")
for res in results:
    print(f"  {res['r']:>4}  {res['dx']:>10.6f}  {res['L1_h']:>12.6e}  {res['L2_h']:>12.6e}  {res['L1_u']:>12.6e}  {res['L2_u']:>12.6e}")

print("\n  ✓ Terminé.\n")