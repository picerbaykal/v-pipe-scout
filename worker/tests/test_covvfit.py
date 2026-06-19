import sys
from pathlib import Path

# Add the project root directory to the Python path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd

from worker.tasks import run_covvfit


#-- helpers ---------------------------------------------------------

class FakeRedis:
    """Capture all .set() calls in-memory."""
    def __init__(self):
        self.store = {}
    def set(self, key, value, ex=None):
        self.store[key] = value

#-- test ------------------------------------------------------------

def test_run_covvfit_minimal(monkeypatch):

    fake_redis = FakeRedis()
    monkeypatch.setattr("worker.tasks.redis_client", fake_redis)
    monkeypatch.setattr(
        "worker.tasks.run_covvfit_inference",
        lambda **kwargs: {"figure_png": "fake_base64_png"},
    )

    # location_data is {location: counts_df}; matrix_df is a DataFrame.
    # Pass real DataFrames directly — the task's isinstance() guards skip
    # deserialization when inputs aren't pickled strings.
    location_data = {"TestCity": pd.DataFrame({"m": [1, 2]})}
    matrix_df = pd.DataFrame([[0.1, 0.9], [0.5, 0.5]])

    # Bound Celery tasks inject `self` automatically when called directly,
    # so we do NOT pass it ourselves.
    result = run_covvfit(location_data, matrix_df)

    assert result == {"figure_png": "fake_base64_png"}