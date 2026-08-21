#!/usr/bin/env python3
"""
LAPIS concurrency saturation test.

Fires N concurrent /sample/aggregated queries at increasing concurrency
levels and measures latency (p50/p95/max) and error rate at each level.
Finds the "knee" — the concurrency where latency starts degrading — so we
can pick MAX_CONCURRENT_CONNECTIONS / BATCH_CONCURRENCY from data instead
of guessing.

Run inside the worker container so it uses the same network path:
    docker compose exec -T worker conda run -n v-pipe-scout-worker \
        python3 /app_shared/lapis_saturation_test.py

Or locally if you can reach LAPIS directly:
    python3 lapis_saturation_test.py
"""
import asyncio
import time
import statistics
import aiohttp

LAPIS_URL = "https://lapis.wasap.genspectrum.org/covid/sample/aggregated"

# a realistic co-occurrence query: two bracketed positions + a date filter
# mirrors what the completeness sweep actually sends
def _make_params(pos_a=22900, pos_b=22920, date="2026-05-01"):
    return {
        "fields": f"[{pos_a}],[{pos_b}]",
        "samplingDate": date,
    }

# concurrency levels to sweep
LEVELS = [1, 2, 4, 8, 16, 24, 32, 48, 64]
# queries per level (enough to see steady-state, not just warmup)
QUERIES_PER_LEVEL = 60
TIMEOUT_S = 30


async def _one_query(session, params):
    t0 = time.monotonic()
    try:
        async with session.get(LAPIS_URL, params=params) as resp:
            await resp.read()
            dt = time.monotonic() - t0
            return (dt, resp.status)
    except asyncio.TimeoutError:
        return (time.monotonic() - t0, "TIMEOUT")
    except Exception as e:
        return (time.monotonic() - t0, f"ERR:{type(e).__name__}")


async def _run_level(concurrency, n_queries):
    """Run n_queries with a semaphore capping concurrency."""
    sem = asyncio.Semaphore(concurrency)
    timeout = aiohttp.ClientTimeout(total=TIMEOUT_S)
    connector = aiohttp.TCPConnector(limit=concurrency, limit_per_host=concurrency)

    async with aiohttp.ClientSession(timeout=timeout, connector=connector) as session:
        async def _guarded(i):
            async with sem:
                # vary positions/date a bit so we don't hit a response cache
                return await _one_query(session, _make_params(
                    pos_a=22900 + (i % 20),
                    pos_b=22920 + (i % 20),
                    date=f"2026-{'%02d' % (3 + i % 4)}-01",
                ))
        t0 = time.monotonic()
        results = await asyncio.gather(*[_guarded(i) for i in range(n_queries)])
        wall = time.monotonic() - t0

    latencies = [r[0] for r in results]
    statuses = [r[1] for r in results]
    ok = [l for l, s in zip(latencies, statuses) if s == 200]
    errors = [s for s in statuses if s != 200]

    return {
        "concurrency": concurrency,
        "n": n_queries,
        "wall_s": wall,
        "throughput_qps": n_queries / wall if wall > 0 else 0,
        "p50": statistics.median(ok) if ok else None,
        "p95": (statistics.quantiles(ok, n=20)[18] if len(ok) >= 20 else max(ok)) if ok else None,
        "max": max(ok) if ok else None,
        "n_ok": len(ok),
        "n_err": len(errors),
        "errors": errors[:5],
    }


async def main():
    print(f"LAPIS saturation test → {LAPIS_URL}")
    print(f"{QUERIES_PER_LEVEL} queries per level, timeout {TIMEOUT_S}s\n")
    print(f"{'conc':>5} {'p50':>8} {'p95':>8} {'max':>8} {'qps':>8} {'ok':>5} {'err':>5}")
    print("-" * 55)

    rows = []
    for level in LEVELS:
        r = await _run_level(level, QUERIES_PER_LEVEL)
        rows.append(r)
        p50 = f"{r['p50']:.2f}s" if r['p50'] else "—"
        p95 = f"{r['p95']:.2f}s" if r['p95'] else "—"
        mx = f"{r['max']:.2f}s" if r['max'] else "—"
        qps = f"{r['throughput_qps']:.1f}"
        print(f"{r['concurrency']:>5} {p50:>8} {p95:>8} {mx:>8} {qps:>8} "
              f"{r['n_ok']:>5} {r['n_err']:>5}")
        if r['errors']:
            print(f"        errors: {r['errors']}")
        # back off between levels so we measure steady-state, not pile-up
        await asyncio.sleep(2)

    # find the knee: last level where p95 stays within 2x of the single-conn p95
    print("\n" + "=" * 55)
    base_p95 = rows[0]["p95"] or rows[0]["p50"]
    if base_p95:
        knee = rows[0]["concurrency"]
        for r in rows:
            if r["n_err"] > 0:
                print(f"⚠ errors appear at concurrency {r['concurrency']}")
                break
            if r["p95"] and r["p95"] <= base_p95 * 2.5:
                knee = r["concurrency"]
        print(f"Baseline p95 (conc=1): {base_p95:.2f}s")
        print(f"Recommended max concurrency (p95 < 2.5x baseline, no errors): {knee}")
        print(f"→ Suggested BATCH_CONCURRENCY per task: {max(2, knee // 2)}")
        print(f"→ Suggested MAX_CONCURRENT_CONNECTIONS (with worker conc 2): {knee}")
        print(f"  (keep total across workers ≤ {knee}: worker_concurrency × BATCH_CONCURRENCY ≤ {knee})")


if __name__ == "__main__":
    asyncio.run(main())