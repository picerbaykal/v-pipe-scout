#!/bin/bash
set -e

echo "=== Skipping pango update — using static pango_summary.json from pango-tree-builder ==="

# 6 processes: one city per process in a 6-city deep scan (2026-10-27);
# set CELERY_CONCURRENCY to change it.
# --max-memory-per-child (kB, 2026-11-04): a process that holds more than this
# after a task is replaced by a fresh one (~250 MB). Without it, a process kept
# 2.3 GB after a deep scan, and the next 6-city scan went over the 6 GB limit
# (OOM, 4 cities lost). Set CELERY_MAX_MEM_PER_CHILD to change it.
echo "=== Starting Celery worker (concurrency ${CELERY_CONCURRENCY:-6}, max ${CELERY_MAX_MEM_PER_CHILD:-600000} kB per process after a task) ==="
exec conda run --no-capture-output -n v-pipe-scout-worker \
    celery -A tasks worker --concurrency "${CELERY_CONCURRENCY:-6}" --prefetch-multiplier 1 \
    --max-memory-per-child "${CELERY_MAX_MEM_PER_CHILD:-600000}" --loglevel=info
