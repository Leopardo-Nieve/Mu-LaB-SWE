"""
mms_sympy.py
Calcule symboliquement la solution manufacturée et les termes sources MMS
pour le solveur Mu-LaB-SWE (SWE 2D stationnaires).

Usage:
    python mms_sympy.py

Sorties:
    - Affichage console des expressions Fortran
    - mms_plot.png : contour plots de h, u, v, S_h, F_x, F_y
"""

import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# Symboles
# ============================================================
x, y   = sp.symbols('x y', real=True)
Lxs    = sp.Symbol('domainX', positive=True)
Lys    = sp.Symbol('domainY', positive=True)
gs     = sp.Symbol('gacl',    positive=True)
nus    = sp.Symbol('nu',      positive=True)

# ============================================================
# Solution manufacturée
#
# Forme choisie :
#   h(x,y) = 0.185 + 0.010 * sin(pi*x/Lx) * sin(2*pi*y/Ly)
#   u(x,y) = 0.124 + 0.030 * sin(pi*x/Lx) * cos(2*pi*y/Ly)
#   v(x,y) =         0.015 * cos(pi*x/Lx) * sin(2*pi*y/Ly)
#
# Propriétés :
#   - sin(pi*x/Lx) = 0 en x=0 et x=Lx → h=0.185, u=0.124, v=0 aux bords
#   - Périodique en y avec période Ly ✓
#   - h > 0 partout (0.175 < h < 0.195) → stable LBM ✓
#   - Perturbations petites (~16% de h0, ~24% de u0) → stabilité numérique ✓
# ============================================================
h_mms = sp.Rational(185,1000) + sp.Rational(10,1000)  * sp.sin(sp.pi*x/Lxs) * sp.sin(2*sp.pi*y/Lys)
u_mms = sp.Rational(124,1000) + sp.Rational(30,1000)  * sp.sin(sp.pi*x/Lxs) * sp.cos(2*sp.pi*y/Lys)
v_mms =                         sp.Rational(15,1000)  * sp.cos(sp.pi*x/Lxs) * sp.sin(2*sp.pi*y/Lys)

# ============================================================
# Termes sources SWE stationnaires (∂/∂t = 0, zb = 0)
#
# Continuité :  ∂(hu)/∂x + ∂(hv)/∂y = S_h
# QDM x :       ∂(huu)/∂x + ∂(huv)/∂y + g*h*∂h/∂x - ν*∇²(hu) = F_x
# QDM y :       ∂(huv)/∂x + ∂(hvv)/∂y + g*h*∂h/∂y - ν*∇²(hv) = F_y
# ============================================================
hu = h_mms * u_mms
hv = h_mms * v_mms

S_h = sp.simplify(sp.diff(hu, x) + sp.diff(hv, y))

F_x = sp.simplify(
      sp.diff(h_mms * u_mms**2, x)
    + sp.diff(h_mms * u_mms * v_mms, y)
    + gs * h_mms * sp.diff(h_mms, x)
    - nus * (sp.diff(hu, x, 2) + sp.diff(hu, y, 2)))

F_y = sp.simplify(
      sp.diff(h_mms * u_mms * v_mms, x)
    + sp.diff(h_mms * v_mms**2, y)
    + gs * h_mms * sp.diff(h_mms, y)
    - nus * (sp.diff(hv, x, 2) + sp.diff(hv, y, 2)))

# ============================================================
# Export Fortran 90
# ============================================================
px = sp.Symbol('position_x')
py = sp.Symbol('position_y')

def to_f90(expr, name):
    e2 = expr.subs({x: px, y: py})
    return sp.fcode(e2, assign_to=name, standard=90,
                    source_format='free', contract=False)

sep = "! " + "="*60

print(sep)
print("! BLOC 1 — à mettre dans MMS_analytic_solution (dans la boucle i,j)")
print(sep)
print(to_f90(h_mms, 'hAnal(i,j)'))
print(to_f90(u_mms, 'uAnal(i,j)'))
print(to_f90(v_mms, 'vAnal(i,j)'))
print()
print(sep)
print("! BLOC 2 — à mettre dans update_body_force (dans la boucle sur la grille centrée)")
print("! Ajouter aux déclarations du module : double precision, allocatable :: S_h_mms(:,:), F_x_mms(:,:), F_y_mms(:,:)")
print(sep)
print(to_f90(S_h, 'S_h_mms(ix,iy)').replace('position_x', 'dx*(DBLE(ix-1)*0.5d0)').replace('position_y', 'dy*(DBLE(iy-1)*0.5d0)'))
print(to_f90(F_x, 'F_x_mms(ix,iy)').replace('position_x', 'dx*(DBLE(ix-1)*0.5d0)').replace('position_y', 'dy*(DBLE(iy-1)*0.5d0)'))
print(to_f90(F_y, 'F_y_mms(ix,iy)').replace('position_x', 'dx*(DBLE(ix-1)*0.5d0)').replace('position_y', 'dy*(DBLE(iy-1)*0.5d0)'))

# ============================================================
# Vérification numérique rapide
# ============================================================
print()
print(sep)
print("! VÉRIFICATION NUMÉRIQUE")
print(sep)
Lx_val, Ly_val = 4.0, 1.0
g_val, nu_val  = 9.81, 1e-4

sub_num = {Lxs: Lx_val, Lys: Ly_val, gs: g_val, nus: nu_val}
h_fn  = sp.lambdify((x, y), h_mms.subs(sub_num), 'numpy')
u_fn  = sp.lambdify((x, y), u_mms.subs(sub_num), 'numpy')
v_fn  = sp.lambdify((x, y), v_mms.subs(sub_num), 'numpy')
Sh_fn = sp.lambdify((x, y), S_h.subs(sub_num),   'numpy')
Fx_fn = sp.lambdify((x, y), F_x.subs(sub_num),   'numpy')
Fy_fn = sp.lambdify((x, y), F_y.subs(sub_num),   'numpy')

xs = np.linspace(0, Lx_val, 50)
ys = np.linspace(0, Ly_val, 40)
XX, YY = np.meshgrid(xs, ys)

print(f"h  min={h_fn(XX,YY).min():.4f}  max={h_fn(XX,YY).max():.4f}  (attendu ~0.185±0.010)")
print(f"u  min={u_fn(XX,YY).min():.4f}  max={u_fn(XX,YY).max():.4f}  (attendu ~0.124±0.030)")
print(f"v  min={v_fn(XX,YY).min():.4f}  max={v_fn(XX,YY).max():.4f}  (attendu ~0.000±0.015)")
print(f"S_h max abs = {np.abs(Sh_fn(XX,YY)).max():.4e}  (petit si h,u,v cohérents)")
print(f"F_x max abs = {np.abs(Fx_fn(XX,YY)).max():.4e}")
print(f"F_y max abs = {np.abs(Fy_fn(XX,YY)).max():.4e}")

# ============================================================
# Plots
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle("MMS — Solution manufacturée et termes sources\n"
             r"$h = 0.185 + 0.01\sin(\pi x/L_x)\sin(2\pi y/L_y)$, "
             r"$u = 0.124 + 0.03\sin(\pi x/L_x)\cos(2\pi y/L_y)$",
             fontsize=12)

fields = [
    (h_fn(XX,YY),  r'$h_{MMS}$ (m)',       'Blues'),
    (u_fn(XX,YY),  r'$u_{MMS}$ (m/s)',      'RdBu_r'),
    (v_fn(XX,YY),  r'$v_{MMS}$ (m/s)',      'RdBu_r'),
    (Sh_fn(XX,YY), r'$S_h$ (source cont.)', 'RdBu_r'),
    (Fx_fn(XX,YY), r'$F_x$ (source QDM x)', 'RdBu_r'),
    (Fy_fn(XX,YY), r'$F_y$ (source QDM y)', 'RdBu_r'),
]

for ax, (field, title, cmap) in zip(axes.flat, fields):
    cf = ax.contourf(XX, YY, field, 20, cmap=cmap)
    plt.colorbar(cf, ax=ax)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')

plt.tight_layout()
plt.savefig('mms_plot.png', dpi=150, bbox_inches='tight')
print("\nPlot sauvegardé : mms_plot.png")