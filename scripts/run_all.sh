#!/usr/bin/env bash
# Compatibility shim: forwards to the Streptococcus pyogenes wrapper
# so older docs/commands that call run_all.sh still work.
set -euo pipefail
exec "$(dirname "$0")/run_all_spa.sh" "$@"
