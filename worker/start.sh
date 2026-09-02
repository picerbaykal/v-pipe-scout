#!/bin/bash
set -e

echo "=== Updating pango signatures ==="

conda run -n v-pipe-scout-worker python3 -c "
import sys
sys.path.insert(0, '/app_shared')
sys.path.insert(0, '/app')
from api.pango_loader import download_pango_summary, PANGO_SUMMARY_CACHE
from pathlib import Path
Path(PANGO_SUMMARY_CACHE).parent.mkdir(parents=True, exist_ok=True)
result = download_pango_summary(PANGO_SUMMARY_CACHE)
if result['success']:
    print(f'Pango updated: {result[\"new_variants\"]} lineages (+{len(result[\"added\"])} new)')
else:
    print(f'Pango update failed: {result[\"error\"]} — using existing file')
"

echo "=== Starting Celery worker ==="
exec conda run --no-capture-output -n v-pipe-scout-worker \
    celery -A tasks worker --concurrency 2 --loglevel=info
