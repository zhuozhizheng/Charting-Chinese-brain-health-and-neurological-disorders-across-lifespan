#!/bin/bash
# Brain-health report agent web service.
# Usage: ./run.sh [port]   (default 8000)
cd "$(dirname "$0")" || exit 1
PORT="${1:-8000}"
CONDA_SH="${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-py38}"
source "$CONDA_SH" && conda activate "$CONDA_ENV" || exit 1
exec python -m brain_agent.cli serve --host 0.0.0.0 --port "$PORT"
