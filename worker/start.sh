#!/bin/bash
set -e

echo "=== Skipping pango update — using static pango_summary.json from pango-tree-builder ==="

echo "=== Starting Celery worker ==="
exec conda run --no-capture-output -n v-pipe-scout-worker \
    celery -A tasks worker --concurrency 2 --loglevel=info
