"""
plot_mms.py
Calcule symboliquement la solution manufacturée et les termes sources MMS
pour le solveur Mu-LaB-SWE (SWE 1D stationnaires, débit uniforme Q).

Manufactured solution (1D in x, uniform in y):
  zb(x)  = 0.3 * sin(4*pi*x/Lx)                   [bed profile from main.f90]
  h(x)   = (2/3)*sin(2*pi*x/Lx) + 2               [free-surface depth]
  u(x)   = Q / h(x)                                [uniform discharge, v=0]
  v(x)   = 0

Continuity is satisfied exactly with no source term (d(hu)/dx = 0 since hu = Q).

Usage:
    python plot_mms.py

Outputs:
    - Console: Fortran 90 code blocks to copy-paste
    - mms_plot.png: contour plots of h, u, v, zb, F_x, F_y
"""

import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

# ============================================================
# Symbols
# ============================================================
x      = sp.Symbol('x', real=True, positive=True)
Lxs    = sp.Symbol('domainX', positive=True)
gs     = sp.Symbol('gacl',    positive=True)
nus    = sp.Symbol('nu',      positive=True)   # substituted directly in Fortran

# ============================================================
# Bed profile  (from main.f90, line 129)
#   zb(x) = 0.3 * sin(4*pi*x / Lx)
# ============================================================
zb_mms    = sp.Rational(3, 10) * sp.sin(4 * sp.pi * x / Lxs)
dzbdx_mms = sp.diff(zb_mms, x)

# ============================================================
# Manufactured solution — uniform discharge Q
#
#   h(x) = (2/3)*sin(2*pi*x/Lx) + 2        [existing Fortran analytic h]
#   u(x) = Q / h(x)                         [d(hu)/dx = 0 exactly]
#   v(x) = 0
#
#   Q = 2 m²/s  →  u(h=2) = 1 m/s  (avoids degenerate zero-velocity case)
# ============================================================
Q = sp.Integer(2)   # m²/s

h_mms = sp.Rational(2, 3) * sp.sin(2 * sp.pi * x / Lxs) + 2
u_mms = Q / h_mms
v_mms = sp.Integer(0)

# ============================================================
# Verify continuity (must be identically zero — no source needed)
# ============================================================
hu = h_mms * u_mms   # = Q (constant)
hv = h_mms * v_mms   # = 0
S_h = sp.simplify(sp.diff(hu, x) + sp.diff(hv, x))
assert S_h == 0, f"Continuity source term is not zero: {S_h}"

# ============================================================
# Steady SWE momentum source terms  (d/dt = 0)
#
#   QDM x:  d(h·u²)/dx  +  g·h·d(h+zb)/dx  −  ν·d²(hu)/dx²  =  F_x
#   QDM y:  0 (v=0, no y-variation)                            =  F_y = 0
#
# Note: hu = Q = const  →  d²(hu)/dx² = 0  (viscous term vanishes)
# ============================================================
F_x = sp.simplify(
      sp.diff(h_mms * u_mms**2, x)
    + gs * h_mms * sp.diff(h_mms + zb_mms, x)
    - nus * sp.diff(h_mms * u_mms, x, 2))

F_y = sp.Integer(0)

# ============================================================
# Fortran 90 code generation
#
# Formatting rules:
#   • double precision literals: sympy's fcode(standard=90) already emits d0
#   • multiples of pi written as  N*pi  (never as 6.28…d0)
#   • trig → double precision:  sin→dsin, cos→dcos
#   • loop indices from source:
#       MMS_analytic_solution  → i, j
#       update_body_force      → ix, iy
#   • nu left as the symbolic variable 'nu' (module variable)
#
# Key insight: sympy fcode(standard=90) already emits Fortran double precision
# notation (e.g. 6.2831853071795865d0).  We replace those exact decimal strings
# with N*pi AFTER the raw fcode call, then apply trig renaming.
# ============================================================
px = sp.Symbol('position_x')

# Map of exact decimal strings sympy emits → Fortran pi-expression
# (sympy uses repr(float(N*pi)) which gives ~16 significant figures in d0 form)
PI_REPLACEMENTS = {
    '25.132741228718345d0': '8.0d0*pi',
    '18.84955592153876d0':  '6.0d0*pi',
    '12.566370614359173d0': '4.0d0*pi',
    '9.42477796076938d0':   '3.0d0*pi',
    '6.2831853071795865d0': '2.0d0*pi',
    '3.1415926535897932d0': 'pi',
    # Negative multiples (rare but covered)
    '-25.132741228718345d0': '-8.0d0*pi',
    '-18.84955592153876d0':  '-6.0d0*pi',
    '-12.566370614359173d0': '-4.0d0*pi',
    '-9.42477796076938d0':   '-3.0d0*pi',
    '-6.2831853071795865d0': '-2.0d0*pi',
    '-3.1415926535897932d0': '-pi',
}

def to_f90_raw(expr, name):
    """Convert sympy expression to Fortran 90 free-form with standard=90 (emits d0)."""
    e2 = expr.subs(x, px)
    return sp.fcode(e2, assign_to=name, standard=90,
                    source_format='free', contract=False)

def fortranify(code):
    """
    Post-process sympy fcode(standard=90) output:
      1. Replace exact pi-decimal strings with N*pi expressions
      2. Rename trig intrinsics to double-precision versions
      3. Remove spurious 'parameter (pi = ...)' lines sympy emits for F_x
      4. Ensure integer literals in expressions carry d0 where needed
    """
    # Step 1: pi multiples (longest/largest first to avoid partial replacement)
    for dec, label in PI_REPLACEMENTS.items():
        code = code.replace(dec, label)

    # Step 2: double-precision trig intrinsics
    code = re.sub(r'\bsin\(', 'dsin(', code)
    code = re.sub(r'\bcos\(', 'dcos(', code)
    code = re.sub(r'\bsqrt\(', 'dsqrt(', code)

    # Step 3: remove the 'parameter (pi = ...)' statement sympy adds for F_x
    code = re.sub(r'^\s*parameter\s*\(pi\s*=.*?\)\s*\n', '',
                  code, flags=re.IGNORECASE | re.MULTILINE)

    # Step 4: replace ALL bare integer literals in the RHS with explicit d0.
    # Reassemble continuation lines first so the partition on '=' is reliable,
    # then split back to free-form lines of ≤80 characters.
    # Bare integer = word-boundary digits NOT already adjacent to '.' or 'd'.
    # Reassemble continuation lines into single logical lines
    # sympy emits ' &\n      ' (space before & and 6-space indent on next line)
    joined = re.sub(r'\s*&\n\s*', '', code)
    fixed_lines = []
    for line in joined.split('\n'):
        if '=' in line and not line.lstrip().startswith('!'):
            lhs, _, rhs = line.partition('=')
            rhs = re.sub(r'(?<![.\d])(\b\d+\b)(?![.\d])', r'\1.0d0', rhs)
            line = lhs + '=' + rhs
        fixed_lines.append(line)
    # Re-wrap to ≤80 chars with Fortran free-form continuation.
    # Break AFTER an operator ( + - * / , ( ) ) so tokens are never split.
    def wrap_f90(line, width=80, indent='      '):
        if len(line) <= width:
            return line
        out = []
        while len(line) > width:
            # Scan backwards from col width-2 for a safe break after an operator
            cut = -1
            for i in range(width - 2, 0, -1):
                if line[i] in ('+', '-', '*', '/', ',', '(', ')'):
                    cut = i + 1   # break AFTER the operator
                    break
            if cut < 1:
                cut = width - 2
            out.append(line[:cut].rstrip() + ' &')
            line = indent + line[cut:].lstrip()
        out.append(line)
        return '\n'.join(out)
    code = '\n'.join(wrap_f90(l) for l in fixed_lines)

    return code

def bloc_analytic(h_expr, u_expr, v_expr):
    """BLOC 1 — MMS_analytic_solution  (loop variables i, j)."""
    parts = []
    for expr, name in [(h_expr, 'hAnal(i,j)'),
                       (u_expr, 'uAnal(i,j)'),
                       (v_expr, 'vAnal(i,j)')]:
        parts.append(fortranify(to_f90_raw(expr, name)))
    return '\n'.join(parts)

def bloc_body_force(Fx_expr, Fy_expr):
    """BLOC 2 — update_body_force  (loop variables ix, iy)."""
    parts = []
    for expr, name in [(Fx_expr, 'force_x_MMS(ix,iy)'),
                       (Fy_expr, 'force_y_MMS(ix,iy)')]:
        parts.append(fortranify(to_f90_raw(expr, name)))
    return '\n'.join(parts)

sep = "! " + "="*62

print(sep)
print("! BLOC 1 — MMS_analytic_solution  (replace hAnal/uAnal/vAnal lines)")
print("!          Loop variables: i, j")
print(sep)
print(bloc_analytic(h_mms, u_mms, v_mms))

print()
print(sep)
print("! BLOC 2 — update_body_force  (replace force_x/y_MMS assignment)")
print("!          Loop variables: ix, iy")
print("!          position_x = dx*(DBLE(ix-1)*0.5d0)  [already set in loop]")
print("!          Uses module variable 'nu' directly — no nu_MMS needed")
print(sep)
print(bloc_body_force(F_x, F_y))

# ============================================================
# Quick numerical verification
# ============================================================
print()
print(sep)
print("! NUMERICAL VERIFICATION")
print(sep)

Lx_val = 4.0
g_val  = 9.81
nu_val = 1.0e-4   # illustrative only (viscous term is analytically zero here)

sub_num = {Lxs: Lx_val, gs: g_val, nus: nu_val}
zb_fn  = sp.lambdify(x, zb_mms.subs(sub_num),  'numpy')
h_fn   = sp.lambdify(x, h_mms.subs(sub_num),   'numpy')
u_fn   = sp.lambdify(x, u_mms.subs(sub_num),   'numpy')
Fx_fn  = sp.lambdify(x, F_x.subs(sub_num),     'numpy')

xs = np.linspace(0, Lx_val, 400)
Q_check = h_fn(xs) * u_fn(xs)
print(f"Q = h*u   min={Q_check.min():.8f}  max={Q_check.max():.8f}  "
      f"(should be constant = {float(Q):.1f})")
print(f"h         min={h_fn(xs).min():.4f}  max={h_fn(xs).max():.4f}  "
      f"(expected: 2 +/- 2/3  ->  [1.333, 2.667])")
print(f"u = Q/h   min={u_fn(xs).min():.4f}  max={u_fn(xs).max():.4f}  "
      f"(varies inversely with h; median ~ 1.0 m/s)")
print(f"F_x       max|.| = {np.abs(Fx_fn(xs)).max():.4e}")

# ============================================================
# Plots
# ============================================================
Ly_val = 1.0
ys = np.linspace(0, Ly_val, 40)
XX, YY = np.meshgrid(xs, ys)

# All fields are 1-D in x — tile into 2-D for contourf
h_2d  = np.tile(h_fn(xs),  (len(ys), 1))
u_2d  = np.tile(u_fn(xs),  (len(ys), 1))
v_2d  = np.zeros_like(h_2d)
zb_2d = np.tile(zb_fn(xs), (len(ys), 1))
Fx_2d = np.tile(Fx_fn(xs), (len(ys), 1))
Fy_2d = np.zeros_like(h_2d)

fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle(
    "MMS: Uniform-discharge solution with sinusoidal bed\n"
    "h = (2/3)*sin(2*pi*x/Lx)+2,  u = Q/h,  v = 0,  Q = 2 m2/s\n"
    "zb = 0.3*sin(4*pi*x/Lx)",
    fontsize=11)

def symplot(ax, field, title, cmap='RdBu_r'):
    """
    Contour plot with colorbar centred on the field median.
    For fields that are analytically uniform (v, F_y), this guarantees
    the colourbar mid-tick sits exactly on the constant value.
    """
    median = float(np.median(field))
    spread = max(abs(float(field.max()) - median),
                 abs(median - float(field.min())))
    # Guard against identically zero or constant fields
    if spread < 1e-10 * (abs(median) + 1e-30):
        spread = 1e-6

    lo = median - spread
    hi = median + spread
    levels = np.linspace(lo, hi, 21)

    cf = ax.contourf(XX, YY, field, levels=levels, cmap=cmap, extend='both')
    cb = plt.colorbar(cf, ax=ax)

    # Ensure median appears as a colourbar tick, removing nearby duplicates
    existing = list(cb.get_ticks())
    tol = spread * 0.05
    cleaned = [t for t in existing if abs(t - median) > tol]
    cb.set_ticks(sorted(cleaned + [median]))
    cb.ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.3g'))

    ax.set_title(title, fontsize=11)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')

plots = [
    (h_2d,  'h_MMS (m)',                        'Blues'),
    (u_2d,  'u_MMS = Q/h (m/s)',                'RdBu_r'),
    (v_2d,  'v_MMS (m/s)  [uniform = 0]',       'RdBu_r'),
    (zb_2d, 'zb bed profile (m)',                'GnBu'),
    (Fx_2d, 'F_x MMS body force (m2/s2)',        'RdBu_r'),
    (Fy_2d, 'F_y MMS body force (m2/s2) [= 0]', 'RdBu_r'),
]

for ax, (field, title, cmap) in zip(axes.flat, plots):
    symplot(ax, field, title, cmap)

plt.tight_layout()
plt.savefig('mms_plot.png', dpi=150, bbox_inches='tight')
print("\nPlot saved: mms_plot.png")
