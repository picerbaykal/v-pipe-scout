"""Implements the Wiseloculus API Queries."""

import logging
import aiohttp
import asyncio
import re
from typing import Optional, List, Tuple, Any, Dict
from datetime import datetime, timedelta

import pandas as pd

from .lapis import Lapis
from .exceptions import APIError
from interface import MutationType

from process.mutations import lapis_mutation_to_pos

# Constants for fallback date range
# When API fails, use the last 3 months instead of an entire year to avoid huge API calls
def get_fallback_date_range() -> Tuple[datetime, datetime]:
    """
    Returns fallback date range of the last 3 months from today.
    This avoids excessively large API calls when the API date range query fails.
    """
    end_date = datetime.now()
    start_date = end_date - timedelta(days=90)  # Approximately 3 months
    return start_date, end_date

FALLBACK_START_DATE, FALLBACK_END_DATE = get_fallback_date_range()

# Connection pool limits to prevent "too many open files" errors
# These limits control the maximum number of concurrent HTTP connections
MAX_CONCURRENT_CONNECTIONS = 100  # Total connections per session
MAX_CONNECTIONS_PER_HOST = 50    # Connections per host

class WiseLoculusLapis(Lapis):
    """Wise-Loculus Instance API"""
    
    @staticmethod
    def _handle_connection_error(error: Exception, context: str = "") -> APIError:
        """
        Convert connection errors to user-friendly APIError messages.
        
        Args:
            error: The exception that occurred
            context: Additional context about where the error occurred
            
        Returns:
            APIError with user-friendly message
        """
        error_msg = str(error)
        
        if "too many open files" in error_msg.lower():
            return APIError(
                f"Too many concurrent connections to the API server. This can happen when querying many locations or a long date range. Try reducing the date range or querying fewer locations at once",
                details=error_msg
            )
        elif "timeout" in error_msg.lower():
            return APIError(
                f"Connection timeout while fetching data from the API server. The server may be busy or unresponsive.",
                details=error_msg
            )
        else:
            prefix = f"Connection error {context}: " if context else "Connection error: "
            return APIError(f"{prefix}{error_msg}", details=error_msg)

    def _mutations_to_and_query(self, mutations: List[str]) -> str:
        """
        Convert a list of mutations to an AND query string for advancedQuery.
        
        Args:
            mutations: List of mutations (e.g., ["23149T", "23224T", "23311T"])
            
        Returns:
            str: AND query string (e.g., "23149T & 23224T & 23311T")
            
        Examples:
            >>> _mutations_to_and_query(["23149T"])
            "23149T"
            >>> _mutations_to_and_query(["23149T", "23224T", "23311T"])
            "23149T & 23224T & 23311T"
            >>> _mutations_to_and_query([])
            ""
        """
        if not mutations:
            return ""
        if len(mutations) == 1:
            return mutations[0]
        return " & ".join(mutations)

    def _transform_query_to_coverage(self, query: str) -> str:
        """
        Transforms a mutation query into a coverage query.
        Replaces each mutation with a check that the position is not N.
        Example: "(S:484K | S:501Y)" -> "(!S:484N | !S:501N)"
        """
        # Regex to capture the mutation parts
        # Group 1: Optional negation (e.g. "!")
        # Group 2: Gene prefix (optional) e.g. "S:", "ORF1a:"
        # Group 3: Position e.g. "501", "23149"
        # Ref base (non-capturing) and Alt base (non-capturing) are ignored
        # Negative lookahead (?!of) prevents matching "3-of" in "[3-of: ...]"
        # Note: No trailing \b because mutations can end with "-" or "." which are not word chars
        pattern = r'(!\s*)?\b([A-Za-z0-9]+:)?(?:[A-Z])?(\d+)[A-Z\-\.](?!of)'
        
        def replace_match(match):
            gene_prefix = match.group(2) or "main:"
            position = match.group(3)
            return f"!{gene_prefix}{position}N"
            
        return re.sub(pattern, replace_match, query)

    async def sample_mutations(
            self, 
            type: MutationType,
            date_range: Tuple[datetime, datetime], 
            locationName: Optional[str] = None,
            min_proportion: float = 0.01,
            nucleotide_mutations: Optional[List[str]] = None,
            amino_acid_mutations: Optional[List[str]] = None,
        ) -> pd.DataFrame:
        """
        Fetches nucleotide mutations for a given date range and optional location.
        Fetches mutations (nucleotide or amino acid) for a given date range and optional location.
        Filters for sequences/reads with particular nucleotide or amino acid mutations, depending on the specified mutation type and provided filters.
        
        Returns a DataFrame with 
        Columns: ['mutation', 'count', 'coverage', 'proportion', 'sequenceName', 'mutationFrom', 'mutationTo', 'position']
        """

        payload = {
            "samplingDateFrom": date_range[0].strftime('%Y-%m-%d'),
            "samplingDateTo": date_range[1].strftime('%Y-%m-%d'),
            "locationName": locationName,
            "minProportion": min_proportion, 
            "orderBy": "proportion",
            "limit": 10000,  # Adjust limit as needed
            "dataFormat": "JSON",
            "downloadAsFile": "false"
        }

        # Add mutation filters if provided
        if nucleotide_mutations:
            payload["nucleotideMutations"] = nucleotide_mutations
        if amino_acid_mutations:
            payload["aminoAcidMutations"] = amino_acid_mutations

        if type == MutationType.AMINO_ACID:
            endpoint = f'{self.server_ip}/sample/aminoAcidMutations'
        elif type == MutationType.NUCLEOTIDE:
            endpoint = f'{self.server_ip}/sample/nucleotideMutations'
        else:
            logging.error(f"Unknown mutation type: {type}")
            return pd.DataFrame()

        try:
            timeout = aiohttp.ClientTimeout(total=30) 
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(
                    endpoint,
                    params=payload,
                    headers={'accept': 'application/json'}
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        df = pd.DataFrame(data['data'])
                        return df
                    else:
                        logging.error(f"Failed to fetch nucleotide mutations: {response.status}")
                        return pd.DataFrame()
        except Exception as e:
            logging.error(f"Error fetching mutations: {e}")
            return pd.DataFrame()

    async def _fetch_mutation_counts_for_date(
            self,
            session: aiohttp.ClientSession,
            locationName: str,
            date_str: str,
            positions: set,
    ) -> List[dict]:
        """
        Fetch mutation counts for a single sampling date via
        /sample/nucleotideMutations, for specific positions only.

        Targeted fetch — only keeps positions from the selected variants'
        pango signatures. No noise filtering needed since we ask for
        specific known positions, not everything circulating.

        Args:
            session: Shared aiohttp session — reused across all dates
                so connection pooling limits apply correctly.
            locationName: Location name e.g. "Lugano (TI)"
            date_str: ISO date string e.g. "2026-01-20"
            positions: Set of integer positions from selected variants'
                pango signatures. Only mutations at these positions
                are kept from the LAPIS response.

        Returns:
            List of dicts with keys: date, pos, cov, var
        """
        params = {
            "locationName": locationName,
            "samplingDateFrom": date_str,
            "samplingDateTo": date_str,
            "minProportion": 0.01,
            "limit": 50000,
        }

        try:
            async with session.get(
                    f'{self.server_ip}/sample/nucleotideMutations',
                    params=params,
                    headers={'accept': 'application/json'}
            ) as response:
                if response.status != 200:
                    error_text = await response.text()
                    raise APIError(
                        f"nucleotideMutations failed for {date_str}: "
                        f"status {response.status}",
                        status_code=response.status,
                        details=error_text,
                        payload=params
                    )
                data = await response.json()
        except APIError:
            raise
        except Exception as e:
            raise APIError(
                f"Error fetching mutations for {date_str}: {str(e)}",
                details=str(e)
            )

        rows = []
        for entry in data.get("data", []):
            # API returns either a 'mutation' field ("C241T") or separate fields
            mutation_str = entry.get("mutation") or ""
            if not mutation_str:
                # reconstruct from separate fields
                mfrom = entry.get("mutationFrom", "")
                mto = entry.get("mutationTo", "")
                position = entry.get("position")
                if position and mto:
                    mutation_str = f"{mfrom}{position}{mto}"
            pos = lapis_mutation_to_pos(mutation_str)
            if not pos:
                continue
            m = re.match(r"^(\d+)", pos)
            if not m or int(m.group(1)) not in positions:
                continue
            coverage = entry.get("coverage", 0)
            if coverage == 0:
                continue
            rows.append({
                "date": date_str,
                "pos": pos,
                "cov": int(coverage),
                "var": int(entry.get("count", 0)),
            })
        return rows

    async def _fetch_mutations_per_date(
            self,
            locationName: str,
            date_range: Tuple[datetime, datetime],
            positions: set,
    ) -> List[dict]:
        """
        Fetch mutation counts across all sampling dates, concurrently.

        First fetches real sampling dates for the location/range (wastewater
        is sampled ~2x/week, not daily), then queries each date in parallel
        using a shared connection-pooled session.

        Args:
            locationName: Location name e.g. "Lugano (TI)"
            date_range: Tuple of (start_date, end_date)
            positions: Set of integer positions from selected variants'
                pango signatures — passed through to each date's fetch.

        Returns:
            List of row dicts {date, pos, cov, var} across all dates.
            Dates that fail are logged and skipped, not fatal.
        """
        dates = await self._get_sampling_dates(locationName, date_range)
        if not dates:
            logging.warning(
                f"No sampling dates found for {locationName} "
                f"{date_range[0].strftime('%Y-%m-%d')} → "
                f"{date_range[1].strftime('%Y-%m-%d')}"
            )
            return []

        connector = aiohttp.TCPConnector(
            limit=MAX_CONCURRENT_CONNECTIONS,
            limit_per_host=MAX_CONNECTIONS_PER_HOST
        )
        timeout = aiohttp.ClientTimeout(total=60)

        async with aiohttp.ClientSession(
                timeout=timeout,
                connector=connector
        ) as session:
            tasks = [
                self._fetch_mutation_counts_for_date(
                    session, locationName, date_str, positions
                )
                for date_str in dates
            ]
            results = await asyncio.gather(*tasks, return_exceptions=True)

        all_rows: List[dict] = []
        failed_dates = []
        for date_str, result in zip(dates, results):
            if isinstance(result, Exception):
                logging.error(
                    f"Failed to fetch mutations for {date_str}: {result}"
                )
                failed_dates.append(date_str)
                continue
            all_rows.extend(result)

        if failed_dates:
            logging.warning(
                f"Mutation fetch failed for {len(failed_dates)}/{len(dates)} "
                f"sampling dates: {failed_dates}"
            )

        logging.info(
            f"Fetched {len(all_rows)} mutation×date rows across "
            f"{len(dates) - len(failed_dates)}/{len(dates)} sampling dates"
        )
        return all_rows

    @staticmethod
    def _rows_to_tallymut(
            rows: List[dict],
            locationName: str,
    ) -> pd.DataFrame:
        """
        Assemble fetched mutation rows into a tallymut-format DataFrame.

        Args:
            rows: List of dicts {date, pos, cov, var} from
                _fetch_mutations_per_date.
            locationName: Location name — added as a column since
                LolliPop expects it in the tallymut.

        Returns:
            pd.DataFrame with tallymut columns:
                date | location | pos | base | cov | var | frac
                + placeholder columns for schema compatibility:
                sample | batch | reads | proto | location_code | gene
        """
        if not rows:
            return pd.DataFrame()

        df = pd.DataFrame(rows)
        df["date"] = pd.to_datetime(df["date"])
        df["location"] = locationName
        df["frac"] = df["var"] / df["cov"]
        df["base"] = df["pos"].str[-1]

        # Placeholder columns required for LolliPop's DataPreprocesser
        # schema — not used in deconvolution computation itself
        df["sample"] = "lapis"
        df["batch"] = "lapis"
        df["reads"] = df["cov"]
        df["proto"] = "lapis"
        df["location_code"] = "0"
        df["gene"] = "genome"

        return df

    async def get_tallymut(
            self,
            locationName: str,
            date_range: Tuple[datetime, datetime],
            variants: List[str],
            pango_loader,
            cowwid_variants=None, #fallback for reconstructed nodes
            reference_positions: set = None,
    ) -> pd.DataFrame:
        """
        Build a tallymut-compatible DataFrame from LAPIS mutation counts,
        for LolliPop deconvolution input.

        Fetches real wastewater mutation counts for the positions defined
        by the selected variants' pango signatures — targeted fetch, no
        noise filtering needed.

        Args:
            locationName: Location name e.g. "Lugano (TI)"
            date_range: Tuple of (start_date, end_date)
            variants: List of pango lineage names selected for the panel.
                Their signatures define which positions to fetch.
            pango_loader: PangoLoader instance — used to get each
                variant's signature mutations.

        Returns:
            pd.DataFrame in tallymut format, ready for LolliPop.
            Empty DataFrame if no data found for this location/range.

        Raises:
            APIError: if LAPIS requests fail.
        """
        # Collect all unique positions from selected variants' signatures
        positions: set = set()
        # Use reference_positions if provided (all cowwid positions)
        # otherwise fall back to selected variants only
        if reference_positions:
            positions = reference_positions
        else:
            for variant in variants:
                if variant in pango_loader._reconstructed_signatures and cowwid_variants and variant in cowwid_variants:
                    signature = cowwid_variants[variant]
                else:
                    signature = pango_loader.get_signature(variant)
                for mut in signature:
                    m = re.match(r"^(\d+)", mut)
                    if m:
                        positions.add(int(m.group(1)))

        if not positions:
            logging.warning(
                f"No positions found for variants {variants} — "
                "check pango_loader has valid signatures."
            )
            return pd.DataFrame()

        logging.info(
            f"Fetching tallymut: {locationName} "
            f"{date_range[0].strftime('%Y-%m-%d')} → "
            f"{date_range[1].strftime('%Y-%m-%d')} | "
            f"{len(variants)} variants | {len(positions)} positions"
        )

        rows = await self._fetch_mutations_per_date(
            locationName, date_range, positions
        )

        df = self._rows_to_tallymut(rows, locationName)

        logging.info(
            f"Built tallymut: {len(df)} rows, "
            f"{df['pos'].nunique() if not df.empty else 0} unique positions, "
            f"{df['date'].nunique() if not df.empty else 0} dates"
        )
        return df

    # ── Co-occurrence fetching ─────────────────────────────────────────────

    async def _fetch_cooccurrence_for_date(
            self,
            session: aiohttp.ClientSession,
            locationName: str,
            date_str: str,
            positions: List[int],
    ) -> List[dict]:
        """
        Fetch read-level co-occurrence at target positions for a single date.

        Positions should be pre-batched to fit within read length (via
        split_positions_by_distance) — otherwise combinations will be
        dominated by rows with N at some position.

        Args:
            session: Shared aiohttp session for connection pooling.
            locationName: e.g. "Lugano (TI)"
            date_str: ISO date string, e.g. "2025-11-09"
            positions: List of positions to query together (typically 2-10,
                       all within max_position_distance_bp).

        Returns:
            List of dicts, one per unique base combination:
                {"date": date_str, "[241]": "T", "[297]": "G", "count": 15000}
        """
        fields = [f"[{p}]" for p in positions]
        params = {
            "locationName": locationName,
            "samplingDateFrom": date_str,
            "samplingDateTo": date_str,
            "fields": ",".join(fields),
        }
        try:
            async with session.get(
                    f"{self.server_ip}/sample/aggregated",
                    params=params,
                    headers={"accept": "application/json"},
            ) as response:
                if response.status != 200:
                    error_text = await response.text()
                    raise APIError(
                        f"cooccurrence failed for {date_str} @ positions={positions}: "
                        f"status {response.status}",
                        status_code=response.status,
                        details=error_text,
                        payload=params,
                    )
                data = await response.json()
        except APIError:
            raise
        except Exception as e:
            raise APIError(
                f"Error fetching cooccurrence for {date_str}: {e}",
                details=str(e),
            )

        rows = []
        for entry in data.get("data", []):
            row = {"date": date_str, "count": int(entry.get("count", 0))}
            for p in positions:
                row[f"[{p}]"] = entry.get(f"[{p}]", "N")
            rows.append(row)
        return rows

    async def get_cooccurrence(
            self,
            locationName: str,
            date_range: Tuple[datetime, datetime],
            positions: List[int],
            dates: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Fetch read-level co-occurrence at target positions across all sampling
        dates in the given range.

        Runs one LAPIS query per sampling date, concurrently. Assumes positions
        have already been batched to fit within read length — call
        split_positions_by_distance first, then call this once per batch.

        Args:
            locationName: Location name (e.g. "Lugano (TI)")
            date_range: (start, end) datetime tuple.
            positions: List of positions in this batch (2-10 recommended,
                       all within max_position_distance_bp).

        Returns:
            DataFrame with columns: date, count, [pos1], [pos2], ...
            One row per (date, unique base combination).
        """
        if dates is None:
            dates = await self._get_sampling_dates(locationName, date_range)
        if not dates:
            logging.warning(
                f"No sampling dates for {locationName} "
                f"{date_range[0].date()} → {date_range[1].date()}"
            )
            return pd.DataFrame()

        connector = aiohttp.TCPConnector(
            limit=MAX_CONCURRENT_CONNECTIONS,
            limit_per_host=MAX_CONNECTIONS_PER_HOST,
        )
        timeout = aiohttp.ClientTimeout(total=120)

        async with aiohttp.ClientSession(timeout=timeout, connector=connector) as session:
            tasks = [
                self._fetch_cooccurrence_for_date(session, locationName, d, positions)
                for d in dates
            ]
            results = await asyncio.gather(*tasks, return_exceptions=True)

        all_rows: List[dict] = []
        failed = []
        for d, res in zip(dates, results):
            if isinstance(res, Exception):
                logging.error(f"cooccurrence fetch failed for {d}: {res}")
                failed.append(d)
                continue
            all_rows.extend(res)

        if failed:
            logging.warning(
                f"cooccurrence failed for {len(failed)}/{len(dates)} dates: {failed[:5]}..."
            )

        logging.info(
            f"Fetched cooccurrence: {len(all_rows)} rows across "
            f"{len(dates) - len(failed)}/{len(dates)} dates | "
            f"{locationName} | positions={positions}"
        )
        return pd.DataFrame(all_rows)

    async def get_queries_over_time(
            self,
            locationName: str,
            date_range: Tuple[datetime, datetime],
            queries: List[Dict[str, str]],
            date_granularity: str = "week",
    ) -> dict:
        """
        Fetch mutation frequencies over time using the queriesOverTime endpoint.

        Args:
            locationName: Location filter e.g. "Lugano (TI)"
            date_range: (start, end) datetime tuple
            queries: List of dicts with keys:
                - countQuery: advanced query string e.g. "main:23018T"
                - coverageQuery: denominator query e.g. "!main:23018N"
                - displayLabel: label for the result (optional)
            date_granularity: "day", "week", or "month"

        Returns:
            Raw response dict with keys:
                queries: list of display labels
                dateRanges: list of {dateFrom, dateTo}
                data: list[query_index][date_range_index] = {count, coverage}
        """
        start, end = date_range

        # build weekly date ranges from start to end
        date_ranges = []
        if date_granularity == "week":
            current = start
            while current <= end:
                week_end = min(current + timedelta(days=6), end)
                date_ranges.append({
                    "dateFrom": current.strftime("%Y-%m-%d"),
                    "dateTo": week_end.strftime("%Y-%m-%d"),
                })
                current = week_end + timedelta(days=1)
        elif date_granularity == "month":
            from calendar import monthrange
            current = start.replace(day=1)
            while current <= end:
                last_day = monthrange(current.year, current.month)[1]
                month_end = min(current.replace(day=last_day), end)
                date_ranges.append({
                    "dateFrom": current.strftime("%Y-%m-%d"),
                    "dateTo": month_end.strftime("%Y-%m-%d"),
                })
                # move to first day of next month
                if current.month == 12:
                    current = current.replace(year=current.year + 1, month=1, day=1)
                else:
                    current = current.replace(month=current.month + 1, day=1)
        else:  # day
            current = start
            while current <= end:
                date_ranges.append({
                    "dateFrom": current.strftime("%Y-%m-%d"),
                    "dateTo": current.strftime("%Y-%m-%d"),
                })
                current += timedelta(days=1)

        payload = {
            "filters": {"locationName": locationName},
            "dateField": "samplingDate",
            "queries": queries,
            "dateRanges": date_ranges,
        }

        async with aiohttp.ClientSession() as session:
            async with session.post(
                    f"{self.server_ip}/component/queriesOverTime",
                    json=payload,
                    headers={"accept": "application/json"},
            ) as response:
                result = await response.json()
                if "error" in result:
                    raise RuntimeError(result["error"])
                return result.get("data", {})


    async def get_date_range(self) -> Tuple[Optional[datetime], Optional[datetime]]:
        """
        Fetches all available sampling dates and returns the earliest and latest dates.
        
        Returns:
            Tuple[Optional[datetime], Optional[datetime]]: (earliest_date, latest_date) or (None, None) if no data
        """
        try:
            timeout = aiohttp.ClientTimeout(total=30)  # 30 second timeout for this query
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(
                    f'{self.server_ip}/sample/aggregated',
                    params={'fields': 'samplingDate'},
                    headers={'accept': 'application/json'}
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        sample_data = data.get('data', [])
                        
                        if not sample_data:
                            logging.warning("No sampling date data available")
                            return None, None
                        
                        # Extract all dates and convert to datetime objects
                        dates = []
                        for entry in sample_data:
                            if 'samplingDate' in entry:
                                try:
                                    date_obj = datetime.strptime(entry['samplingDate'], '%Y-%m-%d')
                                    dates.append(date_obj)
                                except ValueError as e:
                                    logging.warning(f"Invalid date format: {entry['samplingDate']}: {e}")
                        
                        if not dates:
                            logging.warning("No valid sampling dates found")
                            return None, None
                        
                        earliest_date = min(dates)
                        latest_date = max(dates)
                        
                        logging.info(f"Date range: {earliest_date.strftime('%Y-%m-%d')} to {latest_date.strftime('%Y-%m-%d')}")
                        return earliest_date, latest_date
                        
                    else:
                        logging.error(f"Failed to fetch sampling dates: {response.status}")
                        logging.error(await response.text())
                        return None, None
                        
        except Exception as e:
            logging.error(f"Error fetching date range: {e}")
            return None, None

    def get_cached_date_range(self, cache_key: str = "default") -> Tuple[datetime, datetime]:
        """
        Get the date range with caching to avoid repeated API calls.
        Returns pandas Timestamps for compatibility with Streamlit date inputs.
        
        Args:
            cache_key: Unique key for this cache (allows multiple cached ranges)
            
        Returns:
            Tuple[datetime, datetime]: Start and end dates as pandas Timestamps
        """
        import streamlit as st
        import asyncio
        import pandas as pd
        
        # Create a unique session state key
        session_key = f"wiseloculus_date_range_{cache_key}"
        
        # Check if we already have cached date range
        if session_key in st.session_state:
            cached_range = st.session_state[session_key]
            logging.debug(f"Using cached date range for {cache_key}: {cached_range}")
            return cached_range
        
        # Fetch new date range
        try:
            earliest, latest = asyncio.run(self.get_date_range())
            
            if earliest and latest:
                # Convert to pandas Timestamps for Streamlit compatibility
                date_range = (pd.to_datetime(earliest), pd.to_datetime(latest))
                logging.info(f"Fetched date range for {cache_key}: {date_range[0].strftime('%Y-%m-%d')} to {date_range[1].strftime('%Y-%m-%d')}")
            else:
                # Fallback to default dates
                date_range = (pd.to_datetime(FALLBACK_START_DATE), pd.to_datetime(FALLBACK_END_DATE))
                logging.warning(f"API date range not available for {cache_key}, using defaults: {date_range}")
                
        except Exception as e:
            # Fallback to default dates
            date_range = (pd.to_datetime(FALLBACK_START_DATE), pd.to_datetime(FALLBACK_END_DATE))
            logging.warning(f"Error fetching date range for {cache_key}: {e}, using defaults")
        
        # Cache the result
        st.session_state[session_key] = date_range
        return date_range

    def get_cached_date_range_with_bounds(self, cache_key: str = "default") -> Tuple[datetime, datetime, datetime, datetime]:
        """
        Get the date range with bounds for enforcing min/max in date inputs.
        
        Args:
            cache_key: Unique key for this cache
            
        Returns:
            Tuple[datetime, datetime, datetime, datetime]: (start_date, end_date, min_date, max_date)
        """
        start_date, end_date = self.get_cached_date_range(cache_key)
        
        # Use the same dates as bounds to enforce API limits
        # Add a small buffer for edge cases (1 day on each side)
        import pandas as pd
        buffer = pd.Timedelta(days=1)
        min_date = start_date - buffer
        max_date = end_date + buffer
        
        return start_date, end_date, min_date, max_date
    

    async def _component_mutations_over_time(
            self,
            endpoint: str,
            mutation_type_name: str,
            mutations: List[str], 
            date_ranges: List[Tuple[datetime, datetime]],
            locationName: str
        ) -> dict[str, Any]:
        """
        Helper method for fetching mutations over time data from component endpoints.
        
        Args:
            endpoint: The API endpoint name (e.g., "aminoAcidMutationsOverTime")
            mutation_type_name: Display name for logging (e.g., "amino acid")
            mutations: List of mutations
            date_ranges: List of date range tuples
            locationName: Location name to filter by
            
        Returns:
            Dict containing the API response with mutations, dateRanges, and data matrix
        """
        payload = {
            "filters": {
                "locationName": locationName
            },
            "includeMutations": mutations,
            "dateRanges": [
                {
                    "dateFrom": date_range[0].strftime('%Y-%m-%d'),
                    "dateTo": date_range[1].strftime('%Y-%m-%d')
                }
                for date_range in date_ranges
            ],
            "dateField": "samplingDate"
        }

        logging.debug(f"Fetching {mutation_type_name} mutations over time with payload: {payload}")
        
        try:
            timeout = aiohttp.ClientTimeout(total=30)  # 30 second timeout (consistent with other API calls)
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.post(
                    f'{self.server_ip}/component/{endpoint}',
                    headers={
                        'accept': 'application/json',
                        'Content-Type': 'application/json'
                    },
                    json=payload
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    elif response.status == 500:
                        # Log the failed query for debugging
                        logging.error(f"Internal Server Error (500) for {mutation_type_name} mutations over time")
                        logging.error(f"Failed POST query: {self.server_ip}/component/{endpoint}")
                        logging.error(f"Payload: {payload}")
                        error_text = await response.text()
                        logging.error(f"Server response: {error_text}")
                        
                        # Raise custom APIError for better frontend handling
                        raise APIError(
                            f"Internal Server Error: The backend API server is experiencing issues. This is not an application error.",
                            status_code=500,
                            details=error_text,
                            payload=payload
                        )
                    else:
                        logging.error(f"Failed to fetch {mutation_type_name} mutations over time.")
                        logging.error(f"Status code: {response.status}")
                        error_text = await response.text()
                        logging.error(error_text)
                        raise APIError(
                            f"API request failed with status {response.status}",
                            status_code=response.status,
                            details=error_text,
                            payload=payload
                        )
        except APIError:
            # Re-raise our custom APIError
            raise
        except Exception as e:
            logging.error(f"Connection error fetching {mutation_type_name} mutations over time: {e}")
            raise APIError(
                f"Connection error: {str(e)}",
                details=str(e),
                payload=payload
            )

    async def component_aminoAcidMutationsOverTime(
            self, 
            mutations: List[str], 
            date_ranges: List[Tuple[datetime, datetime]],
            locationName: str
        ) -> dict[str, Any]:
        """
        Fetches amino acid mutations over time for a given location and specific date ranges.
        Returns counts and coverage for each mutation and date range.
        
        Args:
            mutations: List of amino acid mutations in format ["S:N501Y", "N:N8N"]
            date_ranges: List of date range tuples [(start_date, end_date), ...]
            locationName: Location name to filter by
            
        Returns:
            Dict containing the API response with mutations, dateRanges, and data matrix
        """
        return await self._component_mutations_over_time(
            endpoint="aminoAcidMutationsOverTime",
            mutation_type_name="amino acid",
            mutations=mutations,
            date_ranges=date_ranges,
            locationName=locationName
        )

    async def component_nucleotideMutationsOverTime(
            self, 
            mutations: List[str], 
            date_ranges: List[Tuple[datetime, datetime]],
            locationName: str
        ) -> dict[str, Any]:
        """
        Fetches nucleotide mutations over time for a given location and specific date ranges.
        Returns counts and coverage for each mutation and date range.
        
        Args:
            mutations: List of nucleotide mutations in format ["A5341C", "C34G"]
            date_ranges: List of date range tuples [(start_date, end_date), ...]
            locationName: Location name to filter by
            
        Returns:
            Dict containing the API response with mutations, dateRanges, and data matrix
        """
        return await self._component_mutations_over_time(
            endpoint="nucleotideMutationsOverTime",
            mutation_type_name="nucleotide",
            mutations=mutations,
            date_ranges=date_ranges,
            locationName=locationName
        )

    def _generate_date_ranges(
            self, 
            date_range: Tuple[datetime, datetime], 
            interval: str = "daily"
        ) -> List[Tuple[datetime, datetime]]:
        """
        Generate date ranges based on the specified interval.
        
        Args:
            date_range: Tuple of (start_date, end_date)
            interval: "daily", "weekly", or "monthly"
            
        Returns:
            List of date range tuples
        """
        start_date, end_date = date_range
        date_ranges = []
        
        if interval == "daily":
            current_date = start_date
            while current_date <= end_date:
                date_ranges.append((current_date, current_date))
                current_date += pd.Timedelta(days=1)
                
        elif interval == "weekly":
            current_date = start_date
            while current_date <= end_date:
                week_end = min(current_date + pd.Timedelta(days=6), end_date)
                date_ranges.append((current_date, week_end))
                current_date = week_end + pd.Timedelta(days=1)
                
        elif interval == "monthly":
            current_date = start_date
            while current_date <= end_date:
                # Get the last day of the current month
                if current_date.month == 12:
                    next_month = current_date.replace(year=current_date.year + 1, month=1, day=1)
                else:
                    next_month = current_date.replace(month=current_date.month + 1, day=1)
                month_end = min(next_month - pd.Timedelta(days=1), end_date)
                date_ranges.append((current_date, month_end))
                current_date = next_month
                
        else:
            raise ValueError(f"Unsupported interval: {interval}. Use 'daily', 'weekly', or 'monthly'")
            
        return date_ranges

    async def mutations_over_time(
            self, 
            mutations: List[str], 
            mutation_type: MutationType, 
            date_range: Tuple[datetime, datetime], 
            locationName: str,
            interval: str = "daily"
        ) -> pd.DataFrame:
        """
        Fetches mutation counts, coverage, and frequency using component endpoints for specified time intervals.

        Args:
            mutations (List[str]): List of mutations to fetch data for.
            mutation_type (MutationType): Type of mutations (NUCLEOTIDE or AMINO_ACID).
            date_range (Tuple[datetime, datetime]): Tuple containing start and end dates for the data range.
            locationName (str): Location name to filter by.
            interval (str): Time interval - "daily" (default), "weekly", or "monthly".

        Returns:
            pd.DataFrame: A MultiIndex DataFrame with mutation and samplingDate as the index, 
                         and count, coverage, and frequency as columns.
        """
        try:
            # Generate date ranges based on the specified interval
            date_ranges = self._generate_date_ranges(date_range, interval)
            
            # Choose the appropriate component endpoint based on mutation type
            if mutation_type == MutationType.AMINO_ACID:
                api_data = await self.component_aminoAcidMutationsOverTime(mutations, date_ranges, locationName)
            elif mutation_type == MutationType.NUCLEOTIDE:
                api_data = await self.component_nucleotideMutationsOverTime(mutations, date_ranges, locationName)
            else:
                raise ValueError(f"Unsupported mutation type: {mutation_type}")

            # Parse the API response (no need to check for "error" key anymore as we raise APIError)
            records = []
            api_data_content = api_data.get("data", {})
            api_mutations = api_data_content.get("mutations", [])
            api_date_ranges = api_data_content.get("dateRanges", [])
            data_matrix = api_data_content.get("data", [])

            # Debug logging
            logging.debug(f"mutations_over_time: Received {len(api_mutations)} mutations, {len(api_date_ranges)} date ranges")

            # Process the data matrix
            for i, mutation in enumerate(api_mutations):
                for j, date_range_info in enumerate(api_date_ranges):
                    if i < len(data_matrix) and j < len(data_matrix[i]):
                        mutation_data = data_matrix[i][j]
                        
                        # Extract count and coverage from the API response
                        count = mutation_data.get("count", 0)
                        coverage = mutation_data.get("coverage", 0)
                        
                        # Calculate frequency from count and coverage
                        frequency = count / coverage if coverage > 0 else pd.NA
                        
                        # For interval-based data, use the start date as the samplingDate
                        # or use the midpoint for better representation
                        start_date = pd.to_datetime(date_range_info["dateFrom"])
                        end_date = pd.to_datetime(date_range_info["dateTo"])
                        
                        if interval == "daily":
                            samplingDate = start_date.strftime('%Y-%m-%d')
                        else:
                            # Use midpoint for weekly/monthly intervals
                            midpoint = start_date + (end_date - start_date) / 2
                            samplingDate = midpoint.strftime('%Y-%m-%d')
                        
                        # Add all records, even those with coverage = 0
                        # This allows us to see when data is missing vs when mutations don't exist
                        records.append({
                            "mutation": mutation,
                            "samplingDate": samplingDate,
                            "count": count,
                            "coverage": coverage,
                            "frequency": frequency
                        })

            logging.debug(f"Created {len(records)} records from component API data")

            # Create DataFrame from records
            df = pd.DataFrame(records)

            # Check for and handle duplicates before setting MultiIndex
            if not df.empty and "mutation" in df.columns and "samplingDate" in df.columns:
                # Check for duplicates
                duplicates = df.duplicated(subset=['mutation', 'samplingDate'], keep=False)
                if duplicates.any():
                    logging.warning(f"Found {duplicates.sum()} duplicate mutation-date combinations, removing duplicates")
                    # Keep first occurrence of each duplicate, preferring non-zero values
                    df = df.drop_duplicates(subset=['mutation', 'samplingDate'], keep='first')
                    logging.debug(f"After deduplication: {len(df)} records remain")
                
                df.set_index(["mutation", "samplingDate"], inplace=True)
            else:
                # Create empty DataFrame with proper MultiIndex structure
                df = pd.DataFrame(columns=["count", "coverage", "frequency"])
                df.index = pd.MultiIndex.from_tuples([], names=["mutation", "samplingDate"])
            return df

        except APIError as api_error:
            # Log the API error and return empty DataFrame
            logging.error(f"APIError encountered in mutations_over_time: {api_error}")
            df = pd.DataFrame(columns=["count", "coverage", "frequency"])
            df.index = pd.MultiIndex.from_tuples([], names=["mutation", "samplingDate"])
            return df
        except Exception as e:
            logging.error(f"Error in mutations_over_time: {e}")
            # Return empty DataFrame with proper MultiIndex structure
            df = pd.DataFrame(columns=["count", "coverage", "frequency"])
            df.index = pd.MultiIndex.from_tuples([], names=["mutation", "samplingDate"])
            return df

    async def coocurrences_over_time(
            self,
            date_range: Tuple[datetime, datetime],
            locationName: str,
            mutations: Optional[List[str]] = None,
            advanced_query: Optional[str] = None,
            interval: str = "daily"
        ) -> pd.DataFrame:
        """
        Fetch proportion data for a SET of mutations (AND filter) or an advanced query over time.
        
        Args:
            date_range: Tuple of (start_date, end_date)
            locationName: Location name to filter by
            mutations: List of nucleotide mutations to filter by (AND condition). Used if advanced_query is None.
            advanced_query: Raw advanced query string. If provided, overrides mutations.
            interval: "daily", "weekly", or "monthly"
            
        Returns:
            DataFrame with columns: samplingDate, count, coverage, frequency 
        """
        try:
            # Determine query string and mutations for coverage
            if advanced_query:
                query_str = advanced_query
                # Transform query to coverage query (replace mutations with !posN)
                coverage_query_str = self._transform_query_to_coverage(advanced_query)
            elif mutations:
                query_str = self._mutations_to_and_query(mutations)
                # For simple mutations list, we can use the same transform logic on the AND query
                # Or use the old method. Let's use the new transform logic for consistency.
                coverage_query_str = self._transform_query_to_coverage(query_str)
            else:
                logging.warning("No mutations or advanced query provided")
                return pd.DataFrame(columns=['samplingDate', 'count', 'coverage', 'frequency'])

            # Step 1: Query for reads matching query and intersection coverage simultaneously
            date_ranges = self._generate_date_ranges(date_range, interval)
            
            filtered_lookup = {}  # date -> count of reads matching query
            coverage_lookup = {}  # date -> reads covering all positions
            
            # Configure TCPConnector with connection limits to prevent "too many open files"
            # This prevents exhausting file descriptors when querying many locations/dates
            connector = aiohttp.TCPConnector(
                limit=MAX_CONCURRENT_CONNECTIONS,
                limit_per_host=MAX_CONNECTIONS_PER_HOST
            )
            
            async with aiohttp.ClientSession(
                timeout=aiohttp.ClientTimeout(total=60),
                connector=connector
            ) as session:
                tasks = []
                
                for date_start, date_end in date_ranges:
                    # Task 1: Filtered count (reads matching query)
                    filtered_payload = {
                        "locationName": locationName,
                        "samplingDateFrom": date_start.strftime('%Y-%m-%d'),
                        "samplingDateTo": date_end.strftime('%Y-%m-%d'),
                        "advancedQuery": query_str,
                        "fields": ["samplingDate"]
                    }
                    
                    filtered_task = session.post(
                        f'{self.server_ip}/sample/aggregated',
                        headers={'accept': 'application/json', 'Content-Type': 'application/json'},
                        json=filtered_payload
                    )
                    
                    # Task 2: Intersection coverage (reads covering all positions)
                    # Require non-N calls at all positions: !posN & !posN ...
                    intersection_payload = {
                        "locationName": locationName,
                        "samplingDateFrom": date_start.strftime('%Y-%m-%d'),
                        "samplingDateTo": date_end.strftime('%Y-%m-%d'),
                        "advancedQuery": coverage_query_str,
                        "fields": ["samplingDate"]
                    }
                    
                    intersection_task = session.post(
                        f'{self.server_ip}/sample/aggregated',
                        headers={'accept': 'application/json', 'Content-Type': 'application/json'},
                        json=intersection_payload
                    )
                    
                    tasks.append((filtered_task, intersection_task))
                
                # Execute all queries in parallel
                all_results = await asyncio.gather(*[task for pair in tasks for task in pair], return_exceptions=True)
                
                # Process results in pairs (filtered, intersection)
                for idx in range(0, len(all_results), 2):
                    filtered_resp = all_results[idx]
                    intersection_resp = all_results[idx + 1]
                    
                    # Process filtered response
                    if isinstance(filtered_resp, Exception):
                        raise self._handle_connection_error(filtered_resp, "fetching filtered data")
                    
                    if filtered_resp.status == 200:
                        filtered_data = await filtered_resp.json()
                        filtered_items = filtered_data.get('data', [])
                        for item in filtered_items:
                            date_str = item.get('samplingDate')
                            count = item.get('count', 0)
                            if date_str:
                                filtered_lookup[date_str] = count
                    else:
                        error_text = await filtered_resp.text()
                        logging.error(f"API error for filtered data: status={filtered_resp.status}, body={error_text}")
                        raise APIError(f"API error (filtered data): {filtered_resp.status}", status_code=filtered_resp.status, details=error_text)

                    # Process coverage (intersection) response
                    if isinstance(intersection_resp, Exception):
                        raise self._handle_connection_error(intersection_resp, "fetching intersection coverage")

                    if intersection_resp.status == 200:
                        intersection_data = await intersection_resp.json()
                        intersection_items = intersection_data.get('data', [])
                        for item in intersection_items:
                            date_str = item.get('samplingDate')
                            count = item.get('count', 0)
                            if date_str:
                                coverage_lookup[date_str] = count
                    else:
                        error_text = await intersection_resp.text()
                        logging.error(f"API error for intersection coverage: status={intersection_resp.status}, body={error_text}")
                        raise APIError(f"API error (coverage data): {intersection_resp.status}", status_code=intersection_resp.status, details=error_text)
            
            # Step 2: Combine coverage and filtered count to calculate frequency
            records = []
            for date_str, coverage in coverage_lookup.items():
                filtered_count = filtered_lookup.get(date_str, 0)
                if coverage > 0:
                    frequency = filtered_count / coverage
                    records.append({
                        'samplingDate': pd.to_datetime(date_str),
                        'count': filtered_count,
                        'coverage': coverage,
                        'frequency': frequency
                    })
            
            # Create DataFrame
            if records:
                df = pd.DataFrame(records)
                df = df.sort_values('samplingDate')
                return df
            else:
                return pd.DataFrame(columns=['samplingDate', 'count', 'coverage', 'frequency'])
                
        except APIError:
            raise
        except OSError as e:
            # Handle OS-level connection errors (e.g., "too many open files") with better messages
            raise self._handle_connection_error(e)
        except aiohttp.ClientError as e:
            # Handle aiohttp-specific errors with better messages
            raise self._handle_connection_error(e)
        except Exception as e:
            logging.error(f"Error in coocurrences_over_time: {e}")
            import traceback
            logging.error(traceback.format_exc())
            raise APIError(f"Unexpected error while fetching co-occurrence data: {str(e)}", details=str(e))

    # ── Tallymut fetching (LAPIS-sourced deconvolution input) ────────────────
    #
    # Used by the Abundance & Co-occurrence tab to source LolliPop deconvolution
    # input directly from LAPIS instead of a pre-built tallymut.tsv file.
    # Wastewater is sampled ~2x/week, we fetch real sampling dates
    # first, then query /sample/nucleotideMutations once per date.

    async def _get_sampling_dates(
        self,
        locationName: str,
        date_range: Tuple[datetime, datetime],
    ) -> List[str]:
        """
        Fetch actual sampling dates available for a location and date range.

        Wastewater is sampled ~2x/week, not daily — this avoids generating
        mostly-empty requests for dates with no samples.

        Returns ISO date strings sorted ascending. Raises APIError on failure
        rather than returning an empty list silently — callers should not
        mistake "API failed" for "no samples exist."
        """

        payload = {
            "locationName": locationName,
            "samplingDateFrom": date_range[0].strftime('%Y-%m-%d'),
            "samplingDateTo": date_range[1].strftime('%Y-%m-%d'),
            "fields": ["samplingDate"],
        }
        try:
            timeout = aiohttp.ClientTimeout(total=30)
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.post(
                        f'{self.server_ip}/sample/aggregated',
                        headers={
                            'accept': 'application/json',
                            'Content-Type': 'application/json'
                        },
                        json=payload
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        dates = sorted({
                            row["samplingDate"]
                            for row in data.get("data", [])
                            if row.get("samplingDate")
                        })
                        logging.info(
                            f"Found {len(dates)} sampling dates for "
                            f"{locationName} "
                            f"{payload['samplingDateFrom']} → "
                            f"{payload['samplingDateTo']}"
                        )
                        return dates
                    else:
                        error_text = await response.text()
                        raise APIError(
                            f"API request failed with status {response.status}",
                            status_code=response.status,
                            details=error_text,
                            payload=payload
                        )
        except APIError:
            raise
        except OSError as e:
            raise self._handle_connection_error(e, "fetching sampling dates")
        except aiohttp.ClientError as e:
            raise self._handle_connection_error(e, "fetching sampling dates")
        except Exception as e:
            logging.error(f"Error fetching sampling dates: {e}")
            raise APIError(
                f"Unexpected error fetching sampling dates: {str(e)}",
                details=str(e)
            )

    