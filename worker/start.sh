#!/bin/bash
set -e

echo "=== Skipping pango update — using static pango_summary.json from pango-tree-builder ==="

# 6 processes: one city per process in a 6-city deep scan (2026-10-27);
# set CELERY_CONCURRENCY to change it
echo "=== Starting Celery worker (concurrency ${CELERY_CONCURRENCY:-6}) ==="
exec conda run --no-capture-output -n v-pipe-scout-worker \
    celery -A tasks worker --concurrency "${CELERY_CONCURRENCY:-6}" --prefetch-multiplier 1 --loglevel=info
