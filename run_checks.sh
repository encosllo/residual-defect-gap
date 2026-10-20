#!/bin/bash
# Reproduces logs/*.log, in 4 parallel processes. Usage: ./run_checks.sh
# The longest job is construction.g (about 80 minutes, mostly Case B).
cd "$(dirname "$0")"
run_job() {
  args="${1%%|*}"; script="${1#*|}"; name="${script%.g}"
  echo "$(date '+%F %T') start $name"
  gap -q -o 6g -c "$args" gap/lib.g "gap/$script" < /dev/null 2>&1 | tr -d '\r' > "logs/$name.log"
  echo "$(date '+%F %T') end   $name"
}
export -f run_job
xargs -P 4 -I{} bash -c 'run_job "$@"' _ {} <<'JOBS'
long:=true;;|construction.g
|smallest.g
|affine.g
|blocks.g
|obstruction.g
|bryce.g
|uniform.g
JOBS
