#!/usr/bin/env bash
# T-String e2e lane (music bead 3ez8g.16): the plucked exact-ZOH modal-string
# fixture through the shared lane engine. See music_lane.sh for the contract.
set -euo pipefail
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/music_lane.sh" --fixture string "$@"
