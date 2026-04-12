#!/usr/bin/env bash
# =============================================================================
# convergence_temporal.sh
# Étude de convergence TEMPORELLE (ordre p en dt) pour Mu-LaB-SWE (LBM)
#
# Usage (depuis Git Bash ou terminal MSYS2 ucrt64):
#   bash convergence_temporal.sh 1 2 4 8
#
# Principe:
#   dx est FIXÉ au niveau le plus fin disponible (dx_fixed)
#   Pour chaque facteur r  →  dt = dt_ref / r
#                          →  tau recalculé pour garder nu constante:
#                               nu  = (tau_ref - 0.5) * dx² / (3 * dt_ref)
#                               tau = 0.5 + 3 * nu * dt / dx²
#                                   = 0.5 + (tau_ref - 0.5) * dt / dt_ref
#                                   = 0.5 + (tau_ref - 0.5) / r
#
# ATTENTION: tau diminue avec r → risque tau < 0.51 pour r grands
#
# Prérequis:
#   - Lancer depuis Git Bash ou terminal MSYS2 ucrt64
#   - Python 3 avec pandas, numpy, matplotlib
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Ajouter gfortran au PATH si nécessaire
# ---------------------------------------------------------------------------
if ! command -v gfortran &>/dev/null; then
    export PATH="/c/msys64/ucrt64/bin:$PATH"
fi
if ! command -v gfortran &>/dev/null; then
    echo "ERREUR: gfortran introuvable."
    exit 1
fi

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SRC_DIR="./src"
BIN_DIR="./bin"
RESULTS_DIR="./results"
MAIN_TEMPLATE="${SRC_DIR}/main.f90"
MODULE_SRC="${SRC_DIR}/Mu_LaB_SWE.f90"
BINARY="${BIN_DIR}/run_convergence.exe"
CONV_RESULTS_DIR="${RESULTS_DIR}/convergence_temporal"

# Paramètres de référence
DX_FIXED=0.01250000   # dx fixé au niveau le plus fin (r=8 de l'étude spatiale)
DT_REF=0.02           # dt de référence (r=1)
TAU_REF=1.982         # tau de référence correspondant à dt_ref et dx_fixed
TAU_MIN=0.51          # seuil de stabilité minimum

# nu physique conservée = (tau_ref - 0.5) * dx² / (3 * dt_ref)
# → calculée dans le script pour vérification

EPSILON_STUDY=1.0d-4
ITERA_MAX=500000

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <r1> <r2> [r3] [r4] ..."
    echo "Exemple: $0 1 2 4 8"
    echo "  → Facteur r: dx fixé=${DX_FIXED}m, dt = ${DT_REF}/r"
    exit 1
fi

REFINEMENTS=("$@")
N_LEVELS=${#REFINEMENTS[@]}

echo "============================================================"
echo " Étude de convergence TEMPORELLE — Mu-LaB-SWE"
echo " dx fixé = ${DX_FIXED} m"
echo " Niveaux de raffinement: ${REFINEMENTS[*]}"
echo "============================================================"

mkdir -p "${BIN_DIR}" "${CONV_RESULTS_DIR}"

# Calcul de nu physique (pour vérification)
NU=$(python3 -c "print(f'{(${TAU_REF} - 0.5) * ${DX_FIXED}**2 / (3.0 * ${DT_REF}):.6e}')")
echo " nu physique conservée = ${NU} m²/s"

# ---------------------------------------------------------------------------
# Boucle principale
# ---------------------------------------------------------------------------
declare -a DT_LIST
declare -a CSV_LIST

for r in "${REFINEMENTS[@]}"; do

    # dt = dt_ref / r
    DT=$(python3 -c "print(f'{${DT_REF} / ${r}:.10f}')")
    # tau = 0.5 + (tau_ref - 0.5) / r   pour conserver nu
    TAU=$(python3 -c "print(f'{0.5 + (${TAU_REF} - 0.5) / ${r}:.6f}')")

    echo ""
    echo "------------------------------------------------------------"
    echo " Niveau r = ${r}"
    printf "   dx   = %s m (fixé)\n" "${DX_FIXED}"
    printf "   dt   = %s s\n" "${DT}"
    printf "   tau  = %s\n"   "${TAU}"
    echo "------------------------------------------------------------"

    # Vérification stabilité tau
    TAU_CHECK=$(python3 -c "t=${TAU}; print('WARN' if t<${TAU_MIN} else 'OK')")
    if [ "${TAU_CHECK}" = "WARN" ]; then
        echo "  ⚠ AVERTISSEMENT: tau=${TAU} < ${TAU_MIN} — instabilité probable, skip ce niveau!"
        continue
    fi

    # --- Patch de main.f90 ---
    PATCHED_MAIN="${SRC_DIR}/main_t${r}.f90"
    cp "${MAIN_TEMPLATE}" "${PATCHED_MAIN}"

    sed -i "s/dx = [0-9]*\.[0-9]*/dx = ${DX_FIXED}/g"                       "${PATCHED_MAIN}"
    sed -i "s/dt = [0-9]*\.[0-9]*d\?[0-9]*/dt = ${DT}d0/g"                  "${PATCHED_MAIN}"
    sed -i "s/tau = [0-9]*\.[0-9]*d\?[0-9]*/tau = ${TAU}d0/g"               "${PATCHED_MAIN}"
    sed -i "s/epsilon = [0-9]*\.[0-9]*d[-+]*[0-9]*/epsilon = ${EPSILON_STUDY}/g" "${PATCHED_MAIN}"
    sed -i "s/itera_no = NINT([^)]*)/itera_no = ${ITERA_MAX}/g"             "${PATCHED_MAIN}"
    sed -i "s/itera_no = [0-9]*/itera_no = ${ITERA_MAX}/g"                  "${PATCHED_MAIN}"

    # --- Compilation ---
    echo "  Compilation..."
    gfortran -O2 "${MODULE_SRC}" "${PATCHED_MAIN}" -o "${BINARY}" 2>&1 | sed 's/^/    [gfortran] /'
    echo "  ✓ Compilation OK"

    # --- Exécution ---
    echo "  Exécution du solveur..."
    "${BINARY}" 2>&1 | sed 's/^/    [run] /'
    # "${BINARY}" 2>&1 | tail -5 | sed 's/^/    [run] /'
    echo "  ✓ Run terminé"

    # --- Récupération du CSV ---
    LATEST_CSV=$(ls -t "${RESULTS_DIR}"/*.csv 2>/dev/null | head -1)
    if [ -z "${LATEST_CSV}" ]; then
        echo "  ERREUR: aucun CSV trouvé dans ${RESULTS_DIR}/"
        exit 1
    fi

    DEST_CSV="${CONV_RESULTS_DIR}/t${r}_dt${DT}.csv"
    cp "${LATEST_CSV}" "${DEST_CSV}"
    echo "  ✓ CSV sauvegardé: ${DEST_CSV}"

    DT_LIST+=("${DT}")
    CSV_LIST+=("${DEST_CSV}")

    rm -f "${PATCHED_MAIN}"

done

echo ""
echo "============================================================"
echo " Runs terminés. Lancer convergence_plot.py pour les graphes."
echo " Résultats dans: ${CONV_RESULTS_DIR}/"
echo "============================================================"