from celery import Celery
import os
import json
import time
import redis
import pickle
import base64
import pandas as pd
from deconvolve import devconvolve
import logging

logger = logging.getLogger(__name__)

from celery.signals import worker_process_init


@worker_process_init.connect
def preload_signatures(**kwargs):
    """Pre-warm signature cache in each worker process at startup."""
    from cooc import get_all_lineage_signatures, get_panel_parent_map
    get_all_lineage_signatures()
    get_panel_parent_map()
    logger.info("Worker process: signatures pre-warmed")


# Initialize Celery
app = Celery(
    'tasks',
    broker=os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0'),
    backend=os.environ.get('CELERY_RESULT_BACKEND', 'redis://redis:6379/0')
)

# Initialize Redis client for storing progress updates
redis_client = redis.Redis(
    host=os.environ.get('REDIS_HOST', 'redis'),
    port=int(os.environ.get('REDIS_PORT', 6379)),
    password=os.environ.get('REDIS_PASSWORD', 'defaultpassword123'),
    db=0
)


# EXAMPLE TASK: that progressively updates its status in Redis
@app.task(bind=True)
def long_running_task(self, n_iterations, sleep_time):
    """
    A simple long-running task that simulates work by sleeping.
    Updates progress in Redis to allow the frontend to track status.
    """
    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    # Initialize result container
    result = {
        "iterations_completed": 0,
        "total_iterations": n_iterations,
        "results": []
    }

    # Process each iteration
    for i in range(n_iterations):
        # Simulate work
        time.sleep(sleep_time)

        # Calculate some dummy result
        iteration_result = {
            "iteration": i + 1,
            "timestamp": time.time(),
            "value": (i + 1) * sleep_time
        }

        # Add to results
        result["results"].append(iteration_result)
        result["iterations_completed"] = i + 1

        # Update progress in Redis
        progress_data = {
            "current": i + 1,
            "total": n_iterations,
            "status": f"Processing iteration {i + 1}/{n_iterations}",
            "partial_results": result["results"]
        }

        redis_client.set(
            progress_key,
            json.dumps(progress_data),
            ex=3600  # Expire after 1 hour
        )

    # Task completed
    progress_data = {
        "current": n_iterations,
        "total": n_iterations,
        "status": "Completed",
        "partial_results": result["results"]
    }

    redis_client.set(
        progress_key,
        json.dumps(progress_data),
        ex=3600  # Expire after 1 hour
    )

    return result


@app.task(bind=True)
def run_deconvolve(self, mutation_counts_df, mutation_variant_matrix_df,
                   bootstraps=None, bandwidth=None, regressor=None,
                   regressor_params=None, deconv_params=None, locationName=None):
    """
    A task that runs the deconvolve function with progress tracking.

    Args:
        mutation_counts_df (pd.DataFrame): DataFrame containing mutation counts data (required)
        mutation_variant_matrix_df (pd.DataFrame): DataFrame containing mutation variant matrix data (required)
        bootstraps (int, optional): Number of bootstrap iterations
        bandwidth (int, optional): Bandwidth parameter for kernel
        regressor (str, optional): Regressor type
        regressor_params (dict, optional): Parameters for the regressor
        deconv_params (dict, optional): Parameters for deconvolution
        locationName (str, optional): Name of the location for proper result structuring
    """
    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    # Initialize progress tracking
    progress_data = {
        "current": 0,
        "total": 5,  # We'll track progress in 5 stages
        "status": "Preparing input data",
        "partial_results": None
    }

    redis_client.set(
        progress_key,
        json.dumps(progress_data),
        ex=3600  # Expire after 1 hour
    )

    try:
        # Update progress
        progress_data["current"] = 1
        progress_data[
            "status"] = f"Preparing deconvolution (bootstraps={bootstraps if bootstraps is not None else 'default'})"
        redis_client.set(progress_key, json.dumps(progress_data), ex=3600)

        # Function to update progress
        def update_progress(stage, message):
            progress_data["current"] = stage
            progress_data["status"] = message
            redis_client.set(progress_key, json.dumps(progress_data), ex=3600)

        # Convert serialized DataFrames back to pandas DataFrames if needed
        try:
            # Add debug info about the input types
            update_progress(1.5,
                            f"Input types: mutation_counts_df: {type(mutation_counts_df)}, mutation_variant_matrix_df: {type(mutation_variant_matrix_df)}")

            # Check if inputs are already DataFrames or need to be deserialized
            if isinstance(mutation_counts_df, pd.DataFrame) and isinstance(mutation_variant_matrix_df, pd.DataFrame):
                update_progress(2, "Inputs are already DataFrames, no parsing needed")
            else:
                # Try to deserialize if needed
                try:
                    # If they're base64 encoded pickle strings
                    if isinstance(mutation_counts_df, str):
                        try:
                            mutation_counts_df = pickle.loads(base64.b64decode(mutation_counts_df))
                            update_progress(2,
                                            f"Successfully unpickled counts DataFrame, shape: {mutation_counts_df.shape}")
                        except:
                            update_progress(2, "Failed to unpickle counts DataFrame as base64")

                    if isinstance(mutation_variant_matrix_df, str):
                        try:
                            mutation_variant_matrix_df = pickle.loads(base64.b64decode(mutation_variant_matrix_df))
                            update_progress(2,
                                            f"Successfully unpickled matrix DataFrame, shape: {mutation_variant_matrix_df.shape}")
                        except:
                            update_progress(2, "Failed to unpickled matrix DataFrame as base64")

                except Exception as e:
                    update_progress(2, f"Error parsing DataFrames: {str(e)}")
                    raise ValueError(f"Failed to deserialize DataFrames: {str(e)}")
        except Exception as e:
            update_progress(2, f"Error processing DataFrames: {str(e)}")
            raise ValueError(f"Failed to process DataFrames: {str(e)}")

        # Create kwargs dict with required parameters and optional parameters if provided
        kwargs = {
            'mutation_counts_df': mutation_counts_df,
            'mutation_variant_matrix_df': mutation_variant_matrix_df
        }

        # Add optional parameters only if they're not None (exclude locationName)
        if bootstraps is not None:
            kwargs['bootstraps'] = bootstraps
        if bandwidth is not None:
            kwargs['bandwidth'] = bandwidth
        if regressor is not None:
            kwargs['regressor'] = regressor
        if regressor_params is not None:
            kwargs['regressor_params'] = regressor_params
        if deconv_params is not None:
            kwargs['deconv_params'] = deconv_params

        # Update progress before running deconvolution
        update_progress(3, "Running deconvolution algorithm")

        # Run the deconvolution with only the provided parameters
        deconvolved_data = devconvolve(**kwargs)

        # Update progress after deconvolution is complete
        update_progress(4, "Processing results")

        # If locationName is provided, restructure the result to use the actual location name
        if locationName and isinstance(deconvolved_data, dict) and "location" in deconvolved_data:
            # Replace the generic "location" key with the actual location name
            location_data = deconvolved_data.pop("location")  # Remove and get the data
            deconvolved_data[locationName] = location_data  # Add with proper name
            update_progress(4.5, f"Restructured result for location: {locationName}")

        # Stage 5: Finalize results
        progress_data["current"] = 5
        progress_data["total"] = 5
        progress_data["status"] = "Completed"
        progress_data["partial_results"] = {"summary": "Deconvolution completed successfully"}
        redis_client.set(progress_key, json.dumps(progress_data), ex=3600)

        return deconvolved_data
    except Exception as e:
        # If there's an error, report it
        error_message = str(e)
        progress_data["status"] = f"Error: {error_message}"
        redis_client.set(progress_key, json.dumps(progress_data), ex=3600)
        raise


@app.task(bind=True)
def run_deconvolve_lapis(self, location: str, start_date: str, end_date: str,
                         variants: list, bootstraps: int = 100, bandwidth: int = 10):
    """
    Celery task for LAPIS-sourced deconvolution.
    Used by the Abundance & Co-occurrence tab.

    Args:
        location: Location name e.g. "Lugano (TI)"
        start_date: ISO date string e.g. "2026-01-06"
        end_date: ISO date string e.g. "2026-04-06"
        variants: List of pango lineage names
        bootstraps: Bootstrap iterations (default 100 = Standard)
        bandwidth: Gaussian kernel bandwidth (default 10 = Narrow)
    """
    from datetime import datetime
    from abundance_cooc import run_deconv_lapis

    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    try:
        redis_client.set(progress_key, json.dumps({
            "current": 1, "total": 4,
            "status": f"Fetching mutation data from LAPIS for {location}..."
        }), ex=3600)

        result = run_deconv_lapis(
            location=location,
            start_date=datetime.fromisoformat(start_date),
            end_date=datetime.fromisoformat(end_date),
            variants=variants,
            bootstraps=bootstraps,
            bandwidth=bandwidth,
            progress_callback=lambda step, msg: redis_client.set(
                progress_key,
                json.dumps({"current": step, "total": 4, "status": msg}),
                ex=3600
            )
        )

        redis_client.set(progress_key, json.dumps({
            "current": 4, "total": 4,
            "status": "Deconvolution complete."
        }), ex=3600)

        return result

    except Exception as e:
        redis_client.set(progress_key, json.dumps({
            "current": 0, "total": 4,
            "status": f"Error: {str(e)}"
        }), ex=3600)
        raise


def _reporter(progress_key: str, total: int = 4, remap=None):
    """Progress writer for the page: {current, total, status, frac, t}.
    frac = how far into the current step (0-1); the page turns it into a % and
    a time left. remap(step) -> step lets a task show an inner function's steps
    as one of its own."""
    import time as _time

    def _w(step, msg, frac=None):
        if remap:
            step = remap(step)
        redis_client.set(progress_key, json.dumps({
            "current": step, "total": total, "status": msg,
            "frac": None if frac is None else round(float(frac), 3),
            "t": _time.time()}), ex=3600)
    return _w


@app.task(bind=True)
def run_cooc_completeness_lapis(self, location: str, start_date: str, end_date: str,
                                variants: list):
    """
    Celery task for LAPIS-sourced panel completeness (cooc).
    Used by the Abundance & Co-occurrence tab.

    Args:
        location: Location name e.g. "Lugano (TI)"
        start_date, end_date: ISO date strings
        variants: List of pango lineage names in the panel
    """
    from datetime import datetime
    from cooc import run_cooc_panel_completeness

    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    try:
        redis_client.set(progress_key, json.dumps({
            "current": 1, "total": 4,
            "status": f"Preparing panel positions for {location}..."
        }), ex=3600)

        result = run_cooc_panel_completeness(
            location=location,
            start_date=datetime.fromisoformat(start_date),
            end_date=datetime.fromisoformat(end_date),
            variants=variants,
            progress_callback=_reporter(progress_key),
        )

        redis_client.set(progress_key, json.dumps({
            "current": 4, "total": 4,
            "status": "Panel completeness computed."
        }), ex=3600)

        return result

    except Exception as e:
        redis_client.set(progress_key, json.dumps({
            "current": 0, "total": 4,
            "status": f"Error: {str(e)}"
        }), ex=3600)
        raise


def _scan_and_name(location, variants, unexplained_patterns, day_totals, position_coverage,
                   hotspots=None):
    """The scanner on unexplained patterns, plus designation dates for the UI."""
    import pandas as pd
    from process.scanner import scan_unexplained_patterns
    from cooc import get_all_lineage_signatures, get_panel_parent_map
    patterns_df = (pd.DataFrame(unexplained_patterns) if unexplained_patterns
                   else pd.DataFrame(columns=["date", "count", "confirmed_present"]))
    result = scan_unexplained_patterns(
        unexplained_patterns=patterns_df, panel_variants=variants,
        all_lineage_signatures=get_all_lineage_signatures(),
        panel_parent_map=get_panel_parent_map(), min_read_count=500,
        day_totals=day_totals, position_coverage=position_coverage, hotspots=hotspots)
    try:
        from api.pango_loader import PangoLoader, get_pango_summary_path as _gp
        _raw = PangoLoader(_gp()).get_raw_data()
        for _c in result.get("resolved_clade", []):
            _c["designation"] = _raw.get(_c["node"], {}).get("designationDate", "")
    except Exception:
        pass
    return result


def _data_positions_cached(location, d0, d1, share, cov):
    """cooc.data_positions, shared between tasks for an hour (every city's
    deep scan needs every location's data for the hotspot positions)."""
    from cooc import data_positions
    key = f"cooc:datapos:{location}:{d0.date()}:{d1.date()}:{share}:{cov}"
    try:
        raw = redis_client.get(key)
        if raw:
            return {int(p): set(a) for p, a in json.loads(raw).items()}
    except Exception:
        pass
    out = data_positions(location, d0, d1, share, cov)
    try:
        redis_client.set(key, json.dumps({str(p): sorted(a) for p, a in out.items()}), ex=3600)
    except Exception:
        pass
    return out


@app.task(bind=True)
def run_cooc_deep_scan_lapis(self, location: str, start_date: str, end_date: str,
                             variants: list, reference_locations: list = None):
    """Phase 2 of a run (2026-10-02): read-level scan at today's positions PLUS
    the positions where the data shows a mutation (cooc.data_positions), then
    the scanner on what the panel doesn't explain. Replaces the scanner task fed
    with phase-1 patterns: those only contained listed positions, so a lineage's
    newest mutations and undesignated ones were never seen. The completeness
    graph stays on phase 1. scope.data_positions: false -> today's positions only.
    """
    import sys
    from datetime import datetime
    sys.path.insert(0, "/app_shared")
    from cooc import run_cooc_panel_completeness, data_positions
    from utils.config import get_cooc_setting
    from process.cooc import _check_cfg

    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    _p = _reporter(progress_key, total=5)
    # the read-level step of run_cooc_panel_completeness (its steps 1-3) is this
    # task's step 2; its step 4 (aggregating) is this task's step 3
    _inner = _reporter(progress_key, total=5, remap=lambda s: 2 if s <= 3 else 3)
    try:
        d0, d1 = datetime.fromisoformat(start_date), datetime.fromisoformat(end_date)
        extra, hot = {}, set()
        if get_cooc_setting("scope.data_positions", default=True):
            _p(1, f"Finding positions with signal in {location}...")
            c = _check_cfg()
            _sh, _cv = float(c["absent_freq"]), int(c["min_cov"])
            extra = _data_positions_cached(location, d0, d1, _sh, _cv)
            # error hotspots: positions where ANY location's data shows >= 2
            # different new bases — an error-prone spot comes from the protocol,
            # it only crosses 1 % in some cities. All available locations, not
            # the run's, so the novel list doesn't depend on the selection.
            hot = {p for p, a in extra.items() if len(a) >= 2}
            refs = [l for l in (reference_locations or []) if l != location]
            for i, ref in enumerate(refs):
                _p(1, f"Error hotspots: {ref}...", (i + 1) / (len(refs) + 1))
                try:
                    hot |= {p for p, a in _data_positions_cached(ref, d0, d1, _sh, _cv).items()
                            if len(a) >= 2}
                except Exception as e:
                    logger.warning(f"[deep scan] hotspots from {ref} failed: {e}")
        _p(2, f"Reading reads (+{len(extra)} positions from the data)...")
        res = run_cooc_panel_completeness(location=location, start_date=d0, end_date=d1,
                                          variants=variants, extra_positions=extra,
                                          progress_callback=_inner)
        _p(3, "Classifying unexplained patterns...")
        day_totals = {str(d)[:10]: int(m) + int(u) for d, m, u in zip(
            res.get("dates", []), res.get("matched_counts", []), res.get("unexplained_counts", []))}
        # error hotspots: positions where the data shows >= 2 different new
        # bases (scanner._drop_hotspots); None without data positions
        hot = sorted(hot) if extra else None
        result = _scan_and_name(location, variants, res.get("unexplained_patterns", []),
                                day_totals, res.get("position_coverage"), hotspots=hot)
        result["deep"] = {"data_positions": len(extra),
                          "hotspots": len(hot) if hot is not None else None}
        # the lineages it found get the panel's own check (★ markers present /
        # measurable per day), so their table cells mean the same as the panel's
        nodes_all = list(dict.fromkeys(
            c["node"] for k in ("resolved_clade", "one_day") for c in (result.get(k) or [])
            if c.get("node")))
        _nmax = int(get_cooc_setting("lists.findings_check_max", default=40))
        nodes = nodes_all[:_nmax]
        if len(nodes_all) > _nmax:
            logger.warning(f"[deep scan] {location}: checking {_nmax} of {len(nodes_all)} "
                           "found lineages (lists.findings_check_max)")
        if nodes:
            from cooc import lineages_check
            _p(4, f"Checking ★ markers of {len(nodes)} found lineages...")
            try:
                result["findings_check"] = lineages_check(location, d0, d1, nodes, variants,
                                                          progress_callback=_p, step=4)
            except Exception as e:
                logger.warning(f"[deep scan] findings check failed in {location}: {e}")
        _p(5, "Deep scan complete.")
        return result
    except Exception as e:
        _p(0, f"Error: {str(e)}")
        raise


@app.task(bind=True)
def run_cooc_lineages_check_lapis(self, location: str, start_date: str, end_date: str,
                                  variants: list, panel: list):
    """Phase 3, cross-check (2026-10-02): the ★ marker check for lineages the
    deep scan named in OTHER cities (cooc.lineages_check). Returns {dates,
    per_variant: {v: {markers, per_date}}}."""
    import sys
    from datetime import datetime
    sys.path.insert(0, "/app_shared")
    from cooc import lineages_check
    _p = _reporter(f"task_progress:{self.request.id}", total=2)
    try:
        _p(1, f"Checking {len(variants)} lineages in {location}...")
        r = lineages_check(location, datetime.fromisoformat(start_date),
                           datetime.fromisoformat(end_date), variants, panel,
                           progress_callback=_p, step=1)
        _p(2, "Done.")
        return r
    except Exception as e:
        _p(0, f"Error: {str(e)}")
        raise


@app.task(bind=True)
def run_cooc_variant_check_lapis(self, location: str, start_date: str, end_date: str,
                                 variant: str, panel: list):
    """"Investigate a variant" -> "Check in data": read counts at one lineage's
    ★ markers in one city (cooc.variant_check)."""
    import sys
    from datetime import datetime
    sys.path.insert(0, "/app_shared")
    from cooc import variant_check
    _p = _reporter(f"task_progress:{self.request.id}", total=2)
    try:
        _p(1, f"Checking {variant} in {location}...")
        r = variant_check(location, datetime.fromisoformat(start_date),
                          datetime.fromisoformat(end_date), variant, panel,
                          progress_callback=_p)
        _p(2, "Done.")
        return r
    except Exception as e:
        _p(0, f"Error: {str(e)}")
        raise


@app.task(bind=True)
def run_cooc_scanner_lapis(self, location: str, start_date: str, end_date: str,
                           variants: list, unexplained_patterns: list,
                           day_totals: dict = None, position_coverage: dict = None):
    """
    Celery task for the co-occurrence panel scanner.

    Takes the unexplained patterns from run_cooc_completeness_lapis and
    classifies them into three buckets:
      - missing_from_panel: cowwid tracked variants not in panel
      - emerging_sublineage: pango descendants of panel variants
      - possibly_new: no known lineage explains the pattern

    Args:
        location: Location name e.g. "Lugano (TI)"
        start_date, end_date: ISO date strings
        variants: Currently selected panel variants
        unexplained_patterns: List of {date, count, confirmed_present} dicts
            from run_cooc_completeness_lapis result["unexplained_patterns"]
        day_totals: {date: matched + unexplained reads} from the same result —
            the denominator of the scanner's per-day evidence share
        position_coverage: {date: {position: reads covering it}} from the same
            result — lets the scanner tell, for a finding seen on one day only,
            whether later samples covered its positions
    """
    import sys, re
    import pandas as pd
    sys.path.insert(0, "/app_shared")

    from api.pango_loader import PangoLoader, get_pango_summary_path
    from api.signatures import get_variant_list
    from process.scanner import scan_unexplained_patterns

    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    try:
        redis_client.set(progress_key, json.dumps({
            "current": 1, "total": 3,
            "status": f"Loading signatures for scanner..."
        }), ex=3600)

        # cowwid surveillance variant names — sigs come from pango_summary.json
        # (single source of truth). yaml sigs are partial and inconsistent
        # with the pango centroid sigs used everywhere else.
        cowwid_names = {v.name for v in get_variant_list().variants}

        redis_client.set(progress_key, json.dumps({
            "current": 2, "total": 3,
            "status": "Loading all pango lineage signatures..."
        }), ex=3600)

        from cooc import get_all_lineage_signatures, get_panel_parent_map
        all_sigs = get_all_lineage_signatures()
        panel_parent_map = get_panel_parent_map()

        redis_client.set(progress_key, json.dumps({
            "current": 3, "total": 3,
            "status": "Classifying unexplained patterns..."
        }), ex=3600)

        patterns_df = (
            pd.DataFrame(unexplained_patterns)
            if unexplained_patterns
            else pd.DataFrame(columns=["date", "count", "confirmed_present"])
        )

        result = scan_unexplained_patterns(
            unexplained_patterns=patterns_df,
            panel_variants=variants,
            all_lineage_signatures=all_sigs,
            panel_parent_map=panel_parent_map,
            min_read_count=500,
            day_totals=day_totals,
            position_coverage=position_coverage,
        )
        # fill designation dates on clade findings for the UI
        try:
            from api.pango_loader import PangoLoader, get_pango_summary_path as _gp
            _raw = PangoLoader(_gp()).get_raw_data()
            for _c in result.get("resolved_clade", []):
                _c["designation"] = _raw.get(_c["node"], {}).get("designationDate", "")
        except Exception:
            pass

        redis_client.set(progress_key, json.dumps({
            "current": 3, "total": 3,
            "status": "Scanner complete."
        }), ex=3600)

        return result

    except Exception as e:
        redis_client.set(progress_key, json.dumps({
            "current": 0, "total": 3,
            "status": f"Error: {str(e)}"
        }), ex=3600)
        raise
