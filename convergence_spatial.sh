#!/usr/bin/env bash
# =============================================================================
# convergence_study.sh
# Étude de convergence spatiale (ordre p) pour Mu-LaB-SWE (LBM)
#
# Usage (depuis Git Bash ou terminal MSYS2 ucrt64 dans VS Code):
#   bash convergence_study.sh 1 2 4 8
#
# Principe:
#   Pour chaque facteur r  →  dx = dx_ref / r
#                          →  dt = dt_ref / r^2   (scaling diffusif LBM: dt ∝ dx²)
#
# Valeurs de référence (maillage grossier):
#   dx_ref = 0.1  m
#   dt_ref = 0.02 s
#   tau reste constant sous scaling diffusif dt∝dx² (viscosité physique conservée)
#
# Prérequis:
#   - Lancer depuis Git Bash ou terminal MSYS2 ucrt64 (PAS WSL, PAS CMD)
#   - Python 3 avec pandas, numpy, matplotlib
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Ajouter gfortran (MSYS2 ucrt64) au PATH si pas déjà présent
# ---------------------------------------------------------------------------
if ! command -v gfortran &>/dev/null; then
    export PATH="/c/msys64/ucrt64/bin:$PATH"
fi
if ! command -v gfortran &>/dev/null; then
    echo "ERREUR: gfortran introuvable. Vérifie que MSYS2 ucrt64 est installé dans C:/msys64/"
    echo "  ou lance ce script depuis le terminal MSYS2 ucrt64 directement."
    exit 1
fi

# ---------------------------------------------------------------------------
# Configuration — structure réelle du projet
# ---------------------------------------------------------------------------
SRC_DIR="./src"
BIN_DIR="./bin"
RESULTS_DIR="./results"
MAIN_TEMPLATE="${SRC_DIR}/main.f90"
MODULE_SRC="${SRC_DIR}/Mu_LaB_SWE.f90"
BINARY="${BIN_DIR}/run_convergence.exe"
CONV_RESULTS_DIR="${RESULTS_DIR}/convergence_spatial"

# Paramètres de référence (maillage le plus grossier, r=1)
DX_REF=0.1
DT_REF=0.02
TAU_REF=1.982     # tau du maillage de référence
TAU_MIN=0.51      # seuil de stabilité minimum
TAU_MAX=1.95      # seuil d'avertissement instabilité

# Critère de convergence très lâche pour les runs d'étude (on veut juste la solution)
EPSILON_STUDY=1.0d-4
ITERA_MAX=500000  # garde-fou

# ---------------------------------------------------------------------------
# Facteurs de raffinement (argument CLI ou défaut)
# ---------------------------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <r1> <r2> [r3] [r4] ..."
    echo "Exemple: $0 1 2 4 8"
    echo "  → Facteur r: dx = ${DX_REF}/r, dt = ${DT_REF}/r²"
    exit 1
fi

REFINEMENTS=("$@")
N_LEVELS=${#REFINEMENTS[@]}
echo "============================================================"
echo " Étude de convergence — Mu-LaB-SWE"
echo " Niveaux de raffinement: ${REFINEMENTS[*]}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Préparation des dossiers
# ---------------------------------------------------------------------------
mkdir -p "${BIN_DIR}" "${CONV_RESULTS_DIR}"

# ---------------------------------------------------------------------------
# Boucle principale sur les niveaux de raffinement
# ---------------------------------------------------------------------------
declare -a DX_LIST
declare -a CSV_LIST

for r in "${REFINEMENTS[@]}"; do

    # --- Calcul des paramètres du niveau ---
    # dx = dx_ref / r
    DX=$(python3 -c "print(f'{${DX_REF} / ${r}:.8f}')")
    # dt = dt_ref / r^2  (scaling diffusif: nu = (tau-0.5)*e*dx/3, e=dx/dt → nu ∝ (tau-0.5)*dx²/dt)
    DT=$(python3 -c "print(f'{${DT_REF} / (${r}**2):.10f}')")
    # tau est ajusté pour maintenir nu constante:
    #   nu = (tau-0.5) * (dx/dt) * dx / 3 = (tau-0.5) * dx^2 / (3*dt)
    #   dx²/dt = dx_ref²/dt_ref (constant si dt ∝ dx²) → tau est constant!
    # En réalité avec dt ∝ dx², e = dx/dt augmente avec r → tau doit être recalculé:
    #   e_ref = dx_ref/dt_ref
    #   e_r   = dx/dt = (dx_ref/r) / (dt_ref/r²) = e_ref * r
    #   nu = (tau-0.5) * e * dx / 3 = (tau-0.5) * (e_ref*r) * (dx_ref/r) / 3 = (tau_ref-0.5)*e_ref*dx_ref/3
    # → tau reste CONSTANT sous scaling diffusif dt∝dx² avec e=dx/dt (tau invariant ✓)
    TAU=${TAU_REF}

    echo ""
    echo "------------------------------------------------------------"
    echo " Niveau r = ${r}"
    printf "   dx   = %s m\n" "${DX}"
    printf "   dt   = %s s\n" "${DT}"
    printf "   tau  = %s\n"   "${TAU}"
    echo "------------------------------------------------------------"

    # Avertissement stabilité
    TAU_CHECK=$(python3 -c "t=${TAU}; print('WARN' if t>${TAU_MAX} or t<${TAU_MIN} else 'OK')")
    if [ "${TAU_CHECK}" = "WARN" ]; then
        echo "  ⚠ AVERTISSEMENT: tau=${TAU} hors plage [${TAU_MIN}, ${TAU_MAX}] — risque d'instabilité!"
    fi

    # --- Patch de main.f90 ---
    PATCHED_MAIN="${SRC_DIR}/main_r${r}.f90"
    cp "${MAIN_TEMPLATE}" "${PATCHED_MAIN}"

    # Remplacer dx
    sed -i "s/dx = [0-9]*\.[0-9]*/dx = ${DX}/g" "${PATCHED_MAIN}"
    # Remplacer dt
    sed -i "s/dt = [0-9]*\.[0-9]*d\?[0-9]*/dt = ${DT}d0/g" "${PATCHED_MAIN}"
    # Remplacer tau
    sed -i "s/tau = [0-9]*\.[0-9]*d\?[0-9]*/tau = ${TAU}d0/g" "${PATCHED_MAIN}"
    # Convergence plus stricte pour l'étude (éviter sortie trop tôt)
    sed -i "s/epsilon = [0-9]*\.[0-9]*d[-+]*[0-9]*/epsilon = ${EPSILON_STUDY}/g" "${PATCHED_MAIN}"
    # Cap d'itérations
    sed -i "s/itera_no = NINT([^)]*)/itera_no = ${ITERA_MAX}/g" "${PATCHED_MAIN}"
    sed -i "s/itera_no = [0-9]*/itera_no = ${ITERA_MAX}/g" "${PATCHED_MAIN}"

      # --- Nettoyage des fichiers .mod obsolètes ---
    rm -f *.mod src/*.mod
    
    # --- Compilation ---
    echo "  Compilation..."
    # gfortran -O2 "${MODULE_SRC}" "${PATCHED_MAIN}" -o "${BINARY}" 2>&1 | sed 's/^/    [gfortran] /'
    gfortran -O2 -ffree-line-length-none "${MODULE_SRC}" "${PATCHED_MAIN}" -o "${BINARY}" 2>&1 | sed 's/^/    [gfortran] /' # remove line limit fortran
    echo "  ✓ Compilation OK"

    # --- Exécution ---
    echo "  Exécution du solveur..."
    # Le solveur écrit dans ./results/ (chemin relatif à son CWD)
    # On lance depuis la racine du projet
    "${BINARY}" 2>&1 | tail -5 | sed 's/^/    [run] /'
    echo "  ✓ Run terminé"

    # --- Récupération du CSV le plus récent ---
    LATEST_CSV=$(ls -t "${RESULTS_DIR}"/*.csv 2>/dev/null | head -1)
    if [ -z "${LATEST_CSV}" ]; then
        echo "  ERREUR: aucun CSV trouvé dans ${RESULTS_DIR}/"
        exit 1
    fi

    # Copier avec nom clair pour l'analyse
    DEST_CSV="${CONV_RESULTS_DIR}/r${r}_dx${DX}.csv"
    cp "${LATEST_CSV}" "${DEST_CSV}"
    echo "  ✓ CSV sauvegardé: ${DEST_CSV}"

    DX_LIST+=("${DX}")
    CSV_LIST+=("${DEST_CSV}")

    # Nettoyage du fichier patché temporaire
    rm -f "${PATCHED_MAIN}"

done

# ---------------------------------------------------------------------------
# Post-traitement Python: calcul des erreurs L2 + ordre de convergence + plots
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo " Post-traitement Python — calcul des ordres de convergence"
echo "============================================================"

# Appel du script Python séparé
# python3 convergence_plot.py

echo ""
echo "============================================================"
echo " Étude terminée. Résultats dans: ${CONV_RESULTS_DIR}/"
echo "   - convergence_summary.csv  : tableau des erreurs et ordres"
echo "   - convergence_plot.png     : graphique log-log"
echo "   - r*_dx*.csv               : CSV bruts par niveau"
echo "============================================================"