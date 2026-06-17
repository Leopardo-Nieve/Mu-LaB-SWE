#!/usr/bin/env bash
# =============================================================================
# convergence_all.sh
# Lance l'étude de convergence spatiale ET temporelle, puis les graphes.
#
# Usage:
#   bash convergence_all.sh
#
# Modifier les listes de facteurs ci-dessous selon tes besoins.
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Facteurs de raffinement à utiliser
# ---------------------------------------------------------------------------
SPATIAL_REFINEMENTS="1 2" # "1 2 4 8"    # dx = dx_ref/r,  dt = dt_ref/r²
TEMPORAL_REFINEMENTS="1 2" # "4 8"   # dx fixé,         dt = dt_ref/r

echo ""
echo "########################################################"
echo "#        CONVERGENCE SPATIALE                         #"
echo "########################################################"
bash convergence_spatial.sh ${SPATIAL_REFINEMENTS}

# echo ""
# echo "########################################################"
# echo "#        CONVERGENCE TEMPORELLE                       #"
# echo "########################################################"
# bash convergence_temporal.sh ${TEMPORAL_REFINEMENTS}

# echo ""
# echo "########################################################"
# echo "#        POST-TRAITEMENT PYTHON                       #"
# echo "########################################################"
# python3 convergence_plot.py

echo ""
echo "########################################################"
echo "  Tout terminé."
echo "########################################################"