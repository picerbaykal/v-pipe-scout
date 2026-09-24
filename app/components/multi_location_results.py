"""
Multi-location results component for V-Pipe Scout.

This module provides a tabbed interface for displaying results from multiple locations,
with progress tracking and per-location visualizations.
"""

import streamlit as st
import json
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import Dict, Any, Optional, Optional
import logging

logger = logging.getLogger(__name__)


def render_location_results_tabs(
        location_tasks: Dict[str, str],
        location_results: Dict[str, Any],
        celery_app,
        redis_client
) -> None:
    """
    Render tabbed interface for multi-location results.

    Args:
        location_tasks: Dictionary mapping location names to task IDs
        location_results: Dictionary mapping location names to result data
        celery_app: Celery application instance
        redis_client: Redis client for progress tracking
    """

    if not location_tasks:
        st.info("No analysis tasks have been started yet.")
        return

    # Create tabs for each location
    locationNames = list(location_tasks.keys())
    tabs = st.tabs([f"📍 {loc}" for loc in locationNames])

    for i, (location, tab) in enumerate(zip(locationNames, tabs)):
        with tab:
            task_id = location_tasks[location]

            # Check if we already have results for this location
            if location in location_results:
                render_single_location_result(location, location_results[location])
            else:
                # Check task status
                render_location_progress(location, task_id, celery_app, redis_client)

    # Show combined download section after all individual location tabs
    st.markdown("---")
    st.subheader("📥 Download Combined Results")

    # Check if all locations have results
    all_complete = len(location_results) == len(location_tasks)

    if all_complete:
        render_combined_download_options(location_results)
    else:
        completed = len(location_results)
        total = len(location_tasks)
        st.info(
            f"⏳ Combined download will be available when all locations are processed ({completed}/{total} complete)")


def render_single_location_result(location: str, result_data: Any) -> None:
    """
    Render results for a single location.

    Args:
        location: Location name
        result_data: Analysis results for the location
    """
    st.success(f"Analysis completed for {location}!")

    if not result_data or len(result_data) == 0:
        st.warning("No results data available to visualize.")
        return
    # Extract variants data based on the actual structure we see
    variants_data = None
    if location in result_data:
        variants_data = result_data[location]
    if not variants_data:
        st.error(f"❌ Could not extract variants data from result for {location}")
        st.write("**Available keys in result_data:**",
                 list(result_data.keys()) if isinstance(result_data, dict) else "Not a dictionary")
        return

    # Create visualization
    try:
        fig = create_variant_plot(variants_data, location)
        if fig:
            # unique key per location so frequent reruns (autorefresh) and the
            # tabbed layout don't create duplicate/greyed chart elements
            st.plotly_chart(fig, use_container_width=True,
                            key=f"deconv_plot_{location}")
        else:
            st.warning("⚠️ Could not create plot - no valid time series data found")
    except Exception as e:
        st.error(f"❌ Error creating plot: {str(e)}")
        return


def _shared_variant_legend(color_map: Dict[str, str]) -> None:
    """One horizontal legend for the whole deconv grid (each plot has
    show_legend=False, so this is the only legend)."""
    items = "".join(
        f"<span style='display:inline-flex;align-items:center;margin:0 14px 4px 0;"
        f"white-space:nowrap;font-size:12px;color:#374151;'>"
        f"<span style='display:inline-block;width:18px;height:3px;border-radius:2px;"
        f"background:{c};margin-right:6px;'></span>{v}</span>"
        for v, c in color_map.items())
    st.markdown(
        "<div style='display:flex;flex-wrap:wrap;align-items:center;"
        "padding:2px 2px 8px;'>" + items + "</div>", unsafe_allow_html=True)


def render_location_grid(location_results: Dict[str, Any], location_names,
                         location_tasks: Dict[str, str],
                         celery_app, redis_client) -> None:
    """Deconvolution plots for ALL selected cities as a small-multiples grid.

    - one shared legend on top (per-plot legends off)
    - consistent variant colours across every city
    - equal plot heights
    - dynamic columns: 1 city -> one big plot; 2-4 -> 2 cols; 5-6 -> 3 cols
    Cities still computing show their progress in-slot.
    """
    names = list(location_names)
    if not names:
        st.info("Select at least one location.")
        return

    cmap = build_variant_color_map(location_results, names)
    if cmap:
        _shared_variant_legend(cmap)

    n = len(names)
    # keep plots readable: at most 2 columns (1 city = one big plot). Overview
    # plots are LINES ONLY (bands muddy at small size); the focus view shows bands.
    ncols = 1 if n == 1 else 2
    height = 480 if n == 1 else 340

    def _one(loc):
        st.markdown(f"<div style='font-size:12.5px;font-weight:600;margin-bottom:2px;'>"
                    f"{loc}</div>", unsafe_allow_html=True)
        if loc in location_results:
            rd = location_results[loc]
            vd = rd.get(loc) if isinstance(rd, dict) and loc in rd else rd
            fig = create_variant_plot(vd, loc, color_map=cmap, show_legend=False,
                                      height=height, title="", show_bands=False)
            if fig is not None:
                fig.update_layout(margin=dict(t=10, b=36, l=48, r=12))
                st.plotly_chart(fig, use_container_width=True,
                                key=f"deconv_grid_{loc}")
            else:
                st.caption("No plottable deconvolution data.")
        elif loc in location_tasks:
            render_location_progress(loc, location_tasks[loc], celery_app, redis_client)
        else:
            st.caption("Not started.")

    with st.container():
        for i in range(0, n, ncols):
            cols = st.columns(ncols)
            for j, loc in enumerate(names[i:i + ncols]):
                with cols[j]:
                    _one(loc)

    # ── focus one city at full size (with confidence bands) — a reliable zoom,
    # since the plotly modebar zoom is fiddly on small grid plots ──
    ready = [l for l in names if l in location_results]
    if len(ready) > 1:
        st.markdown("<div style='margin-top:6px'></div>", unsafe_allow_html=True)
        focus = st.selectbox("🔍 Focus one city (full size, with confidence bands)",
                             ["—"] + ready, key="acooc_deconv_focus")
        if focus and focus != "—":
            rd = location_results[focus]
            vd = rd.get(focus) if isinstance(rd, dict) and focus in rd else rd
            fig = create_variant_plot(vd, focus, color_map=cmap, show_legend=True,
                                      height=520, title=f"{focus}", show_bands=True)
            if fig is not None:
                st.plotly_chart(fig, use_container_width=True,
                                key=f"deconv_focus_{focus}")


def render_location_progress(location: str, task_id: str, celery_app, redis_client) -> None:
    """
    Render progress for a location that's still processing.

    Args:
        location: Location name
        task_id: Celery task ID
        celery_app: Celery application instance
        redis_client: Redis client
    """

    # Check if task is completed
    try:
        task = celery_app.AsyncResult(task_id)
        if task.ready():
            try:
                result = task.get()
                # Store result in session state
                st.session_state.location_results[location] = result
                st.success(f"Analysis completed for {location}!")
                st.rerun()
            except Exception as e:
                st.error(f"Error retrieving result for {location}: {str(e)}")
                return
    except Exception as e:
        st.error(f"Error checking task status for {location}: {str(e)}")
        return

    # Show progress
    progress_key = f"task_progress:{task_id}"
    try:
        progress_data = redis_client.get(progress_key)

        if progress_data:
            progress_info = json.loads(progress_data)
            current = progress_info.get('current', 0)
            total = progress_info.get('total', 1)
            status = progress_info.get('status', 'Processing...')

            # Display progress bar
            progress_value = current / total if total > 0 else 0
            st.progress(progress_value)
            st.write(f"Status: {status}")

            if current > 0 and total > 0:
                st.caption(f"Progress: {current}/{total} ({progress_value:.1%})")
        else:
            st.info("Task is running... Progress information will appear shortly.")
    except Exception as e:
        logger.warning(f"Error retrieving progress for {location}: {e}")
        st.info("Task is running... Progress information unavailable.")

    # Add manual check button
    if st.button(f"🔄 Check Status for {location}", key=f"check_{location}"):
        st.rerun()


def _fallback_palette():
    return px.colors.qualitative.Bold if hasattr(px.colors.qualitative, 'Bold') else [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    ]


def build_variant_color_map(location_results: Dict[str, Any],
                            location_names) -> Dict[str, str]:
    """One colour per variant, STABLE across every city, so a shared legend and
    small-multiples grid stay consistent (create_variant_plot otherwise colours
    by per-city index, which would give the same variant different colours in
    different cities). Variants are ordered first-seen across cities;
    'undetermined' is always grey and last."""
    seen = []
    for loc in location_names:
        rd = location_results.get(loc) or {}
        vd = rd.get(loc) if isinstance(rd, dict) and loc in rd else rd
        if isinstance(vd, dict):
            for v in vd.keys():
                if v != "undetermined" and v not in seen:
                    seen.append(v)
    palette = _fallback_palette()
    palette = palette * (len(seen) // len(palette) + 1)
    cmap = {v: palette[i] for i, v in enumerate(seen)}
    cmap["undetermined"] = "#9ca3af"
    return cmap


def create_variant_plot(variants_data: Dict[str, Any], location: str,
                        color_map: Optional[Dict[str, str]] = None,
                        show_legend: bool = True,
                        height: int = 500,
                        title: Optional[str] = None,
                        show_bands: bool = True) -> Optional[go.Figure]:
    """
    Create plotly figure for variant abundance over time.

    Args:
        variants_data: Dictionary containing variant time series data
        location: Location name for the title
        color_map: optional {variant: colour} for consistent colours across cities
        show_legend: draw the per-figure legend (False when a shared legend is used)
        height: figure height in px
        title: figure title override (None -> default per-location title)

    Returns:
        Plotly figure object or None if no valid data
    """
    if not variants_data:
        logger.warning(f"No variants data provided for {location}")
        return None

    fig = go.Figure()

    # Color palette for variants (fallback when no color_map is given)
    colors = _fallback_palette()

    # Ensure we have enough colors
    if len(variants_data) > len(colors):
        # Extend by cycling through the list
        colors = colors * (len(variants_data) // len(colors) + 1)

    # Track the max time range for x-axis limits
    all_dates = []
    traces_added = 0

    # Plot each variant
    for i, (variant_name, variant_data) in enumerate(variants_data.items()):
        logger.info(f"Processing variant: {variant_name}, data type: {type(variant_data)}")

        if not isinstance(variant_data, dict):
            logger.warning(f"Variant {variant_name} data is not a dict: {type(variant_data)}")
            continue

        # Look for time series data with different possible key names
        timeseries = None
        timeseries_key_found = None
        for key in ['timeseriesSummary', 'timeseries', 'time_series', 'summary', 'data', 'points']:
            if key in variant_data:
                timeseries = variant_data[key]
                timeseries_key_found = key
                break

        if not timeseries:
            logger.warning(
                f"No timeseries data found for variant {variant_name}. Available keys: {list(variant_data.keys())}")
            # Try to see if the variant_data itself might be the timeseries
            if isinstance(variant_data, list):
                timeseries = variant_data
                timeseries_key_found = "direct_list"
            else:
                continue

        if not isinstance(timeseries, list):
            logger.warning(
                f"Timeseries for {variant_name} is not a list (found under '{timeseries_key_found}'): {type(timeseries)}")
            continue

        if len(timeseries) == 0:
            logger.warning(f"Empty timeseries for {variant_name}")
            continue

        logger.info(
            f"Found timeseries for {variant_name} under key '{timeseries_key_found}' with {len(timeseries)} points")

        try:
            # Extract dates and proportions
            dates = []
            proportions = []
            lower_bounds = []
            upper_bounds = []

            for j, point in enumerate(timeseries):
                if not isinstance(point, dict):
                    logger.warning(f"Point {j} for {variant_name} is not a dict: {type(point)}")
                    continue

                # Look for date field
                date_val = None
                date_key_found = None
                for date_key in ['date', 'time', 'timestamp', 'Date', 'Time']:
                    if date_key in point:
                        date_val = point[date_key]
                        date_key_found = date_key
                        break

                if date_val is None:
                    logger.warning(
                        f"No date field found in point {j} for {variant_name}. Available keys: {list(point.keys())}")
                    continue

                # Look for proportion/value field
                prop_val = None
                prop_key_found = None
                for prop_key in ['proportion', 'value', 'abundance', 'frequency', 'Proportion', 'Value']:
                    if prop_key in point:
                        prop_val = point[prop_key]
                        prop_key_found = prop_key
                        break

                if prop_val is None:
                    logger.warning(
                        f"No proportion field found in point {j} for {variant_name}. Available keys: {list(point.keys())}")
                    continue

                try:
                    parsed_date = pd.to_datetime(date_val)
                    dates.append(parsed_date)
                    proportions.append(float(prop_val))

                    # Look for confidence intervals
                    lower_val = point.get('proportionLower', point.get('lower', point.get('Lower', prop_val)))
                    upper_val = point.get('proportionUpper', point.get('upper', point.get('Upper', prop_val)))
                    lower_bounds.append(float(lower_val))
                    upper_bounds.append(float(upper_val))

                except (ValueError, TypeError) as e:
                    logger.warning(
                        f"Error parsing data point {j} for {variant_name}: {e} (date='{date_val}', prop='{prop_val}')")
                    continue

            if not dates:
                logger.warning(f"No valid data points found for variant {variant_name}")
                continue

            all_dates.extend(dates)
            if color_map and variant_name in color_map:
                color = color_map[variant_name]
            else:
                color = colors[i % len(colors)]

            # Add line plot for this variant
            fig.add_trace(go.Scatter(
                x=dates,
                y=proportions,
                mode='lines+markers',
                line=dict(color=color, width=2),
                name=variant_name,
                customdata=list(zip(lower_bounds, upper_bounds)),  # Pass confidence intervals as custom data
                hovertemplate=f'<b>{variant_name}</b><br>' +
                              '%{y:.1%} [%{customdata[0]:.1%}, %{customdata[1]:.1%}]<extra></extra>'
            ))

            # Add shaded confidence interval if we have different upper/lower bounds
            # (skipped in the small-multiples grid — the overlapping bands are
            # muddy at small size; the focused single-city view shows them)
            if show_bands and any(l != u for l, u in zip(lower_bounds, upper_bounds)):
                rgba_color = get_rgba_color(color, 0.2)

                fig.add_trace(go.Scatter(
                    x=dates + dates[::-1],  # Forward then backwards
                    y=upper_bounds + lower_bounds[::-1],  # Upper then lower bounds
                    fill='toself',
                    fillcolor=rgba_color,
                    line=dict(color='rgba(0,0,0,0)'),
                    hoverinfo="skip",
                    showlegend=False
                ))

            traces_added += 1
            logger.info(f"Successfully added trace for {variant_name} with {len(dates)} data points")

        except Exception as e:
            logger.error(f"Error processing variant {variant_name}: {e}")
            continue

    if traces_added == 0:
        logger.warning(f"No valid traces added for location {location}")
        return None

    # Update layout
    fig.update_layout(
        title=(title if title is not None
               else f"Variant Proportion Estimates - {location}"),
        xaxis_title="Date",
        yaxis_title="Estimated Proportion",
        yaxis=dict(
            tickformat='.0%',  # Format as percentage
            range=[0, 1]
        ),
        legend_title="Variants",
        showlegend=show_legend,
        height=height,
        margin=dict(t=42, b=40, l=50, r=16),
        template="plotly_white",
        hovermode="x unified"
    )

    logger.info(f"Successfully created plot for {location} with {traces_added} variants")
    return fig


def get_rgba_color(color: str, alpha: float) -> str:
    """
    Convert a color to RGBA format with specified alpha.

    Args:
        color: Color string (hex, rgb, or named color)
        alpha: Alpha value (0-1)

    Returns:
        RGBA color string
    """
    try:
        if color.startswith('#'):
            # Convert hex to rgb
            r = int(color[1:3], 16) / 255
            g = int(color[3:5], 16) / 255
            b = int(color[5:7], 16) / 255
            return f'rgba({r:.3f}, {g:.3f}, {b:.3f}, {alpha})'
        elif color.startswith('rgb'):
            # Extract RGB values from the string
            rgb_values = color.replace('rgb(', '').replace('rgba(', '').replace(')', '').split(',')
            if len(rgb_values) >= 3:
                r = float(rgb_values[0].strip()) / 255 if float(rgb_values[0].strip()) > 1 else float(
                    rgb_values[0].strip())
                g = float(rgb_values[1].strip()) / 255 if float(rgb_values[1].strip()) > 1 else float(
                    rgb_values[1].strip())
                b = float(rgb_values[2].strip()) / 255 if float(rgb_values[2].strip()) > 1 else float(
                    rgb_values[2].strip())
                return f'rgba({r:.3f}, {g:.3f}, {b:.3f}, {alpha})'

        # Fallback for unknown color formats
        return f'rgba(0.5, 0.5, 0.5, {alpha})'
    except Exception:
        # Safe fallback
        return f'rgba(0.5, 0.5, 0.5, {alpha})'


def render_download_options(location: str, variants_data: Dict[str, Any]) -> None:
    """
    Render download buttons for location-specific results.

    Args:
        location: Location name
        variants_data: Variant analysis results
    """

    # Prepare CSV data
    all_variant_data = []

    for variant_name, variant_data in variants_data.items():
        if not isinstance(variant_data, dict):
            continue

        # Look for time series data with different possible key names
        timeseries = None
        for key in ['timeseriesSummary', 'timeseries', 'time_series', 'summary', 'data', 'points']:
            if key in variant_data:
                timeseries = variant_data[key]
                break

        # Fallback: check if variant_data itself is a list
        if not timeseries and isinstance(variant_data, list):
            timeseries = variant_data

        if not timeseries or not isinstance(timeseries, list):
            continue

        for point in timeseries:
            if not isinstance(point, dict):
                continue

            # Extract data with flexible field names
            date_val = None
            for date_key in ['date', 'time', 'timestamp', 'Date', 'Time']:
                if date_key in point:
                    date_val = point[date_key]
                    break

            prop_val = None
            for prop_key in ['proportion', 'value', 'abundance', 'frequency', 'Proportion', 'Value']:
                if prop_key in point:
                    prop_val = point[prop_key]
                    break

            if date_val is None or prop_val is None:
                continue

            all_variant_data.append({
                'location': location,
                'variant': variant_name,
                'date': str(date_val),
                'proportion': prop_val,
                'proportionLower': point.get('proportionLower', point.get('lower', '')),
                'proportionUpper': point.get('proportionUpper', point.get('upper', ''))
            })

    if all_variant_data:
        st.subheader(f"📥 Download {location} Results")
        st.caption(f"Download results for {location} only")

        col1, col2 = st.columns(2)

        with col1:
            csv_data = pd.DataFrame(all_variant_data).to_csv(index=False)
            st.download_button(
                label=f"📄 Download {location} (CSV)",
                data=csv_data,
                file_name=f'deconvolution_results_{location.replace(" ", "_")}.csv',
                mime='text/csv',
                key=f"csv_download_{location}",
                help=f"CSV format results for {location} only"
            )

        with col2:
            json_data = json.dumps({location: variants_data}, indent=2)
            st.download_button(
                label=f"📋 Download {location} (JSON)",
                data=json_data,
                file_name=f'deconvolution_results_{location.replace(" ", "_")}.json',
                mime='application/json',
                key=f"json_download_{location}",
                help=f"JSON format results for {location} only"
            )
    else:
        st.warning(f"No data available for download from {location}")


def render_combined_download_options(location_results: Dict[str, Any]) -> None:
    """
    Render download options for combined multi-location results.

    Args:
        location_results: Dictionary of location -> result data
    """
    if not location_results:
        return

    st.caption("Download all location results in a single long-format table with location column")

    # Prepare combined data for download
    combined_data = []

    for location, result_data in location_results.items():
        # Handle different result data structures
        if isinstance(result_data, dict):
            # Check if result_data has location as a key (prefer actual location name)
            if location in result_data:
                variants_data = result_data[location]
            elif "location" in result_data:
                # Handle generic "location" key structure (legacy)
                variants_data = result_data["location"]
            else:
                # Result data directly contains variants
                variants_data = result_data
        else:
            logger.warning(f"Unexpected result data type for {location}: {type(result_data)}")
            continue

        if not isinstance(variants_data, dict):
            logger.warning(f"Variants data for {location} is not a dict: {type(variants_data)}")
            continue

        # Extract timeseries data for each variant
        for variant_name, variant_data in variants_data.items():
            if not isinstance(variant_data, dict):
                continue

            # Look for time series data with different possible key names
            timeseries = None
            timeseries_key_found = None
            for key in ['timeseriesSummary', 'timeseries', 'time_series', 'summary', 'data', 'points']:
                if key in variant_data:
                    timeseries = variant_data[key]
                    timeseries_key_found = key
                    break

            if not timeseries:
                # Try to see if the variant_data itself might be the timeseries
                if isinstance(variant_data, list):
                    timeseries = variant_data
                    timeseries_key_found = "direct_list"
                else:
                    continue

            if not isinstance(timeseries, list):
                continue

            # Process each time point
            for point in timeseries:
                if not isinstance(point, dict):
                    continue

                # Extract data with flexible field names
                date_val = None
                for date_key in ['date', 'time', 'timestamp', 'Date', 'Time']:
                    if date_key in point:
                        date_val = point[date_key]
                        break

                prop_val = None
                for prop_key in ['proportion', 'value', 'abundance', 'frequency', 'Proportion', 'Value']:
                    if prop_key in point:
                        prop_val = point[prop_key]
                        break

                if date_val is None or prop_val is None:
                    continue

                # Add to combined data with location information
                combined_data.append({
                    'location': location,
                    'variant': variant_name,
                    'date': str(date_val),
                    'proportion': prop_val,
                    'proportionLower': point.get('proportionLower', point.get('lower', '')),
                    'proportionUpper': point.get('proportionUpper', point.get('upper', ''))
                })

    # Show download options
    if combined_data:
        # Download buttons
        col1, col2 = st.columns(2)

        with col1:
            csv_data = pd.DataFrame(combined_data).to_csv(index=False)
            st.download_button(
                label="📄 Download Combined Results (CSV)",
                data=csv_data,
                file_name='deconvolution_results_all_locations.csv',
                mime='text/csv',
                key="combined_csv_download",
                help="Long-format table with location column for all results"
            )

        with col2:
            json_data = json.dumps(location_results, indent=2)
            st.download_button(
                label="📋 Download Combined Results (JSON)",
                data=json_data,
                file_name='deconvolution_results_all_locations.json',
                mime='application/json',
                key="combined_json_download",
                help="Complete nested JSON structure with all results"
            )
    else:
        st.warning("No data available for combined download")