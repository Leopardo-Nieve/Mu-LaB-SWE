"""
convergence_plot.py
Post-traitement des études de convergence spatiale ET temporelle — Mu-LaB-SWE
 
Usage:
    python convergence_plot.py
 
Détecte automatiquement:
    - ./results/convergence_spatial/  → CSV r*_dx*.csv  (convergence spatiale)
    - ./results/convergence_temporal/ → CSV t*_dt*.csv  (convergence temporelle)
 
Sorties dans chaque dossier:
    - convergence_summary.csv  (L1, L2, ordres pour h et u)
    - convergence_plot.png
Et à la racine de results/:
    - convergence_plot_combined.png
"""
 
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
 
COL_H     = "h (m)"
COL_U     = "u (m/s)"
COL_HANAL = "h analytical (m)"
COL_UANAL = "u analytical (m/s)"
 
STUDIES = [
    {
        "dir":      "./results/convergence_spatial",
        "pattern":  re.compile(r"^r\d+_dx([\d.]+)\.csv$"),
        "x_label":  "Δx (m)",
        "x_name":   "dx",
        "title":    "Convergence spatiale",
        "p_expect": 2,
    },
    # {
    #     "dir":      "./results/convergence_temporal",
    #     "pattern":  re.compile(r"^t\d+_dt([\d.]+)\.csv$"),
    #     "x_label":  "Δt (s)",
    #     "x_name":   "dt",
    #     "title":    "Convergence temporelle",
    #     "p_expect": 2,
    # },
]
 
# ---------------------------------------------------------------------------
def compute_orders(x_vals, errors):
    orders = []
    for i in range(len(x_vals) - 1):
        r = np.log(x_vals[i] / x_vals[i+1])
        orders.append(np.log(errors[i] / errors[i+1]) / r if errors[i+1] > 0 else np.nan)
    return orders
 
def analyze(study):
    d, pattern = study["dir"], study["pattern"]
    x_label, x_name = study["x_label"], study["x_name"]
    title, p_expect = study["title"], study["p_expect"]
 
    if not os.path.isdir(d):
        print(f"\n  [SKIP] Dossier introuvable: {d}")
        return None
 
    entries = []
    for fname in os.listdir(d):
        m = pattern.match(fname)
        if m:
            entries.append((float(m.group(1)), os.path.join(d, fname)))
 
    if len(entries) < 2:
        print(f"\n  [SKIP] Moins de 2 fichiers trouvés dans {d}/")
        return None
 
    entries.sort(key=lambda e: -e[0])   # grossier → fin
    x_vals    = np.array([e[0] for e in entries])
    n         = len(entries)
 
    print(f"\n{'='*70}")
    print(f"  {title} ({n} niveaux) :")
    for xv, path in entries:
        print(f"    {x_name}={xv:.8f}  →  {os.path.basename(path)}")
    print(f"{'='*70}")
 
    l1_h, l1_u = [], []
    l2_h, l2_u = [], []
 
    for _, csv_path in entries:
        df = pd.read_csv(csv_path, skipinitialspace=True)
        df.columns = df.columns.str.strip()
 
        # print("columns =", list(df.columns))
        # print(df.head(3).to_string())
 
        # h_num  = df[COL_H].values
        # u_num  = df[COL_U].values
        # h_anal = df[COL_HANAL].values
        # u_anal = df[COL_UANAL].values
 
        h_num  = df["zb (m)"].values
        u_num  = df["h (m)"].values
        h_anal = df["q (m^2/s)"].values
        u_anal = df["h analytical (m)"].values
 
        err_h = np.abs(h_num - h_anal)
        err_u = np.abs(u_num - u_anal)
 
        # print(f"\nDEBUG FILE: {csv_path}")
        # print("u_num[:5]   =", u_num[:5])
        # print("u_anal[:5]  =", u_anal[:5])
        # print("err_u[:5]   =", err_u[:5])
        # print("mean err_u  =", np.mean(err_u))
        # print("max err_u   =", np.max(err_u))
 
        n_pts = len(h_num)
        l1_h.append(np.sum(err_h) / n_pts)
        l1_u.append(np.sum(err_u) / n_pts)
        l2_h.append(np.sqrt(np.sum(err_h**2) / n_pts))
        l2_u.append(np.sqrt(np.sum(err_u**2) / n_pts))
 
    l1_h, l1_u = np.array(l1_h), np.array(l1_u)
    l2_h, l2_u = np.array(l2_h), np.array(l2_u)
 
    ord_l1_h = compute_orders(x_vals, l1_h)
    ord_l1_u = compute_orders(x_vals, l1_u)
    ord_l2_h = compute_orders(x_vals, l2_h)
    ord_l2_u = compute_orders(x_vals, l2_u)
 
    # --- Tableau console ---
    print(f"\n  Norme L1:")
    print(f"  {x_name:>10}  {'L1(h)':>14}  {'p_L1(h)':>9}  {'L1(u)':>14}  {'p_L1(u)':>9}")
    print(f"  {'-'*62}")
    for i in range(n):
        ph = f"{ord_l1_h[i-1]:9.3f}" if i > 0 else "        —"
        pu = f"{ord_l1_u[i-1]:9.3f}" if i > 0 else "        —"
        print(f"  {x_vals[i]:10.6f}  {l1_h[i]:14.6e}  {ph}  {l1_u[i]:14.6e}  {pu}")
 
    print(f"\n  Norme L2:")
    print(f"  {x_name:>10}  {'L2(h)':>14}  {'p_L2(h)':>9}  {'L2(u)':>14}  {'p_L2(u)':>9}")
    print(f"  {'-'*62}")
    for i in range(n):
        ph = f"{ord_l2_h[i-1]:9.3f}" if i > 0 else "        —"
        pu = f"{ord_l2_u[i-1]:9.3f}" if i > 0 else "        —"
        print(f"  {x_vals[i]:10.6f}  {l2_h[i]:14.6e}  {ph}  {l2_u[i]:14.6e}  {pu}")
 
    for label, orders_h, orders_u in [("L1", ord_l1_h, ord_l1_u), ("L2", ord_l2_h, ord_l2_u)]:
        vph = [p for p in orders_h if not np.isnan(p)]
        vpu = [p for p in orders_u if not np.isnan(p)]
        if vph: print(f"\n  Ordre moyen p_{label}(h) = {np.mean(vph):.3f}  (attendu ≈ {p_expect})")
        if vpu: print(f"  Ordre moyen p_{label}(u) = {np.mean(vpu):.3f}")
 
    # --- Sauvegarde CSV résumé ---
    rows = [{
        x_name:      x_vals[i],
        "L1_h":      l1_h[i], "L1_u": l1_u[i],
        "L2_h":      l2_h[i], "L2_u": l2_u[i],
        "order_L1_h": ord_l1_h[i-1] if i > 0 else None,
        "order_L1_u": ord_l1_u[i-1] if i > 0 else None,
        "order_L2_h": ord_l2_h[i-1] if i > 0 else None,
        "order_L2_u": ord_l2_u[i-1] if i > 0 else None,
    } for i in range(n)]
    summary_path = os.path.join(d, "convergence_summary.csv")
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    print(f"\n  Résumé sauvegardé : {summary_path}")
 
   
 
    return dict(
        x_vals=x_vals,
        l1_h=l1_h, l1_u=l1_u,
        l2_h=l2_h, l2_u=l2_u,
        ord_l1_h=ord_l1_h, ord_l1_u=ord_l1_u,
        ord_l2_h=ord_l2_h, ord_l2_u=ord_l2_u,
        x_label=x_label, title=title, outdir=d,
    )
 
# ---------------------------------------------------------------------------
def plot_study(axes_row, res):
    """
    axes_row : liste de 4 axes [L1_h, L1_u, L2_h, L2_u]
    """
    x_vals = res["x_vals"]
    mid    = len(x_vals) // 2
    x_ref  = np.array([x_vals.min() * 0.8, x_vals.max() * 1.2])
 
    configs = [
        (axes_row[0], res["l1_h"], res["ord_l1_h"], "h (profondeur)", "L1", "#2196F3"),
        (axes_row[1], res["l1_u"], res["ord_l1_u"], "u (vitesse)",    "L1", "#F44336"),
        (axes_row[2], res["l2_h"], res["ord_l2_h"], "h (profondeur)", "L2", "#1565C0"),
        (axes_row[3], res["l2_u"], res["ord_l2_u"], "u (vitesse)",    "L2", "#B71C1C"),
    ]
 
    for ax, errors, orders, var_name, norm, color in configs:
        ax.loglog(x_vals, errors, 'o-', color=color, lw=2, ms=7,
                  label=f"Erreur {norm} — {var_name}")
 
        for p_ref, ls, lbl in [(1, '--', 'pente 1'), (2, ':', 'pente 2')]:
            scale = errors[mid] / x_vals[mid]**p_ref
            ax.loglog(x_ref, scale * x_ref**p_ref, ls=ls, color='gray', lw=1.2, label=lbl)
 
        for i in range(len(x_vals) - 1):
            p = orders[i]
            if not np.isnan(p):
                xm = np.sqrt(x_vals[i] * x_vals[i+1])
                ym = np.sqrt(errors[i] * errors[i+1])
                ax.annotate(f"p={p:.2f}", xy=(xm, ym), fontsize=9, color=color,
                            ha='center', va='bottom',
                            bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.75))
 
        ax.set_xlabel(res["x_label"], fontsize=10)
        ax.set_ylabel(f"Erreur {norm}", fontsize=10)
        ax.set_title(f"{res['title']} — {norm}({var_name})", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, which='both', ls=':', alpha=0.5)
 
# ---------------------------------------------------------------------------
results = [analyze(s) for s in STUDIES]
results = [r for r in results if r is not None]
 
if not results:
    print("\nAucune étude à tracer.")
    exit(0)
 
# Plots individuels (1 fichier par étude, 2×2 : L1/L2 × h/u)
for res in results:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"Étude de convergence — {res['title']} (LBM D2Q9)",
                 fontsize=13, fontweight='bold')
    plot_study([axes[0][0], axes[0][1], axes[1][0], axes[1][1]], res)
    plt.tight_layout()
    path = os.path.join(res["outdir"], "convergence_plot.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Plot sauvegardé : {path}")
    plt.close()
 
# Plot combiné si les deux études sont disponibles (4 lignes × 2 colonnes)
if len(results) == 2:
    fig, axes = plt.subplots(4, 2, figsize=(13, 18))
    fig.suptitle("Études de convergence — Mu-LaB-SWE (LBM D2Q9)",
                 fontsize=14, fontweight='bold')
    for row_base, res in enumerate(results):
        # L1 sur la ligne row_base*2, L2 sur row_base*2+1
        plot_study([
            axes[row_base*2][0],   axes[row_base*2][1],
            axes[row_base*2+1][0], axes[row_base*2+1][1],
        ], res)
    plt.tight_layout()
    combined_path = "./results/convergence_plot_combined.png"
    plt.savefig(combined_path, dpi=150, bbox_inches='tight')
    print(f"  Plot combiné sauvegardé : {combined_path}")
    plt.close()
 
print("\n  ✓ Post-traitement terminé.\n")