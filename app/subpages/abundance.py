"""This page allows to compose the list of variants and their respective mutational signatures.
    Steps: 
    1. Select the variants of interest with their respective signature mutations / or load a pre-defined signature mutation set
        1.1 For a variant, search for the signature mutaitons
        1.2 Or load a pre-defined signature mutation set
    2. Build the mutation-variant matrix
    3. Visualize the mutation-variant matrix
    4. Export and download the mutation-variant matrix and var_dates.yaml
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import plotly.graph_objects as go
from pydantic import BaseModel
from typing import List
import logging
import os
import json
import pickle  
import base64 
import asyncio
from datetime import datetime
from celery import Celery
import redis

from interface import MutationType
from api.signatures import get_variant_list, get_variant_names
from api.signatures import Mutation
from api.covspectrum import CovSpectrumLapis
from api.wiseloculus import WiseLoculusLapis
from components.variant_signature_component import render_signature_composer
from state import AbundanceEstimatorState
from api.signatures import Variant as SignatureVariant
from api.signatures import VariantList as SignatureVariantList
from process.variants import create_mutation_variant_matrix
from utils.config import get_wiseloculus_url, get_covspectrum_url
from utils.url_state import create_url_state_manager


# Initialize Celery
celery_app = Celery(
    'tasks',
    broker=os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0'),
    backend=os.environ.get('CELERY_RESULT_BACKEND', 'redis://redis:6379/0')
)

# Initialize Redis client for checking task status
redis_client = redis.Redis(
    host=os.environ.get('REDIS_HOST', 'redis'),
    port=int(os.environ.get('REDIS_PORT', 6379)),
    password=os.environ.get('REDIS_PASSWORD', 'defaultpassword123'),
    db=0
)


class Variant(BaseModel):
    """
    Model for a variant with its signature mutations.
    This is a simplified version of the Variant class from signatures.py.
    """
    name: str  # pangolin name
    signature_mutations: List[str]
    
    @classmethod
    def from_signature_variant(cls, signature_variant: SignatureVariant) -> "Variant":
        """Convert a signature Variant to our simplified Variant."""
        return cls(
            name=signature_variant.name,  # This is already the pangolin name
            signature_mutations=signature_variant.signature_mutations
        )


class VariantList(BaseModel):
    """Model for a simplified list of variants."""
    variants: List[Variant] = []
    
    @classmethod
    def from_signature_variant_list(cls, signature_variant_list: SignatureVariantList) -> "VariantList":
        """Convert a signature VariantList to our simplified VariantList."""
        variant_list = cls()
        for signature_variant in signature_variant_list.variants:
            variant_list.add_variant(Variant.from_signature_variant(signature_variant))
        return variant_list
        
    def add_variant(self, variant: Variant):
        self.variants.append(variant)
        
    def remove_variant(self, variant: Variant):
        self.variants.remove(variant)


@st.cache_data(ttl=3600)  # Cache for 1 hour
def cached_get_variant_list() -> SignatureVariantList:
    """Cached version of get_variant_list to avoid repeated API calls."""
    return get_variant_list()

@st.cache_data(ttl=3600)  # Cache for 1 hour
def cached_get_variant_names() -> List[str]:
    """Cached version of get_variant_names to avoid repeated API calls."""
    return get_variant_names()

def app():
    # Initialize URL state manager for this page
    url_state = create_url_state_manager("abundance")
    
    # ============== INITIALIZATION ==============
    # Initialize all session state variables
    AbundanceEstimatorState.initialize()
    
    # Apply clearing flag for manual inputs if needed
    AbundanceEstimatorState.apply_clear_flag()
    
    # Create a reference to the persistent variant list
    combined_variants = AbundanceEstimatorState.get_combined_variants()
    
    # Check if we need to handle first-time loading of variants
    # If the combined list is empty and we have selected curated names, register them
    if not combined_variants.variants and AbundanceEstimatorState.get_selected_curated_names():
        selected_names = AbundanceEstimatorState.get_selected_curated_names()
        all_curated_variants = cached_get_variant_list().variants
        curated_variant_map = {v.name: v for v in all_curated_variants}
        
        from state import VariantSource
        for name in selected_names:
            if name in curated_variant_map and not AbundanceEstimatorState.is_variant_registered(name):
                variant_to_add = curated_variant_map[name]

                AbundanceEstimatorState.register_variant(
                    name=variant_to_add.name,
                    signature_mutations=variant_to_add.signature_mutations,
                    source=VariantSource.CURATED
                )
        
        # Refresh the combined variants after registering
        combined_variants = AbundanceEstimatorState.get_combined_variants()
    
    # ============== UI HEADER ==============
    # Now start the UI
    st.title("Rapid Variant Abundance Estimation")
    st.subheader("Compose the list of variants and their respective mutational signatures, then estimate their abundance in recent wastewater.")
    st.write("This page allows you to select variants of interest and their respective signature mutations.")
    st.write("You can either select from a curated list, compose a custom variant signature with live queries to CovSpectrum or manually input the mutations.")
    st.write("The selected variants will be used to build a mutation-variant matrix.")
    st.write("Fetch counts and coverage for these mutations.")
    st.write("And finally, estimate the abundance of the variants in the wastewater data.")

    st.markdown("---")

    # ============== CONFIGURATION ==============
    # Setup APIs from centralized config
    cov_spectrum_api = get_covspectrum_url()
    covSpectrum = CovSpectrumLapis(cov_spectrum_api)

    # ============== VARIANT SELECTION ==============
    st.subheader("Variant Selection")
    st.write("Select the variants of interest from either the curated list or compose new signature on the fly from CovSpectrum.")
    
    # combined_variants is already initialized from session state

    # --- Curated Variants Section ---
    st.markdown("#### Curated Variant List")
    # Get the available variant names from the signatures API (cached)
    available_variants = cached_get_variant_names()
    
    # Load variant selection from URL or use current session state
    url_selected_variants = url_state.load_from_url("selected_variants", AbundanceEstimatorState.get_selected_curated_names(), list)
    
    # Filter URL variants to only include those that are still available
    current_selected_curated = AbundanceEstimatorState.get_selected_curated_names()
    available_variant_names = set(available_variants)
    filtered_url_variants = [v for v in url_selected_variants if v in available_variant_names] if url_selected_variants else []
    
    # Use current session state as default, but respect filtered URL state if it exists and differs
    if filtered_url_variants != current_selected_curated and filtered_url_variants:
        # URL has different variants, use filtered URL state
        default_variants = filtered_url_variants
    else:
        # Use current session state, but only include variants that are still available
        default_variants = [v for v in current_selected_curated if v in available_variant_names]
    
    # Create a multi-select box for variants
    selected_curated_variants = st.multiselect(
        "Select known variants of interest – curated by the V-Pipe team",
        options=available_variants,
        default=default_variants,
        help="Select from the list of known variants. The signature mutations of these variants have been curated by the V-Pipe team",
        placeholder="Start typing to search for variants..."
    )
    
    # Save variant selection to URL
    url_state.save_to_url(selected_variants=selected_curated_variants)
    
    # Update the session state if the selection has changed
    if selected_curated_variants != AbundanceEstimatorState.get_selected_curated_names():
        AbundanceEstimatorState.set_selected_curated_names(selected_curated_variants)
        
        # Force immediate registration of selected variants
        all_curated_variants = cached_get_variant_list().variants
        curated_variant_map = {v.name: v for v in all_curated_variants}
        
        # First, unregister any curated variants that are no longer selected
        from state import VariantSource
        registered_variants = AbundanceEstimatorState.get_registered_variants()
        for variant_name, variant_data in list(registered_variants.items()):
            if variant_data['source'] == VariantSource.CURATED and variant_name not in selected_curated_variants:
                AbundanceEstimatorState.unregister_variant(variant_name)
        
        # Then register all newly selected variants
        for name in selected_curated_variants:
            if name in curated_variant_map and not AbundanceEstimatorState.is_variant_registered(name):
                variant_to_add = curated_variant_map[name]
                signature_variant = Variant.from_signature_variant(variant_to_add)
                AbundanceEstimatorState.register_variant(
                    name=signature_variant.name,
                    signature_mutations=signature_variant.signature_mutations,
                    source=VariantSource.CURATED
                )
        
        st.rerun()
    
    # Sync the combined_variants with selected curated variants
    # First, remove any curated variants that are no longer selected
    all_curated_variants = cached_get_variant_list().variants
    curated_names = {v.name for v in all_curated_variants}
    
    # Get custom variant names to avoid removing them
    from state import VariantSource
    custom_variant_names = {
        variant_data['name'] 
        for variant_data in AbundanceEstimatorState.get_registered_variants().values()
        if variant_data['source'] in [VariantSource.CUSTOM_COVSPECTRUM, VariantSource.CUSTOM_MANUAL]
    }
    
    # Create a copy of the list to safely iterate and remove
    variants_to_remove = []
    for v in combined_variants.variants:
        # Only remove curated variants that are no longer selected
        # Make sure NOT to remove custom variants here
        if v.name in curated_names and v.name not in selected_curated_variants and v.name not in custom_variant_names:
            variants_to_remove.append(v)
    
    # Now perform the removal
    for v in variants_to_remove:
        # Remove from the registry first - this is the single source of truth
        AbundanceEstimatorState.unregister_variant(v.name)
        
    # Then add newly selected curated variants if they're not already in the combined list
    if selected_curated_variants:
        # Build a map of name to variant object for quick lookup
        existing_variant_names = {v.name for v in combined_variants.variants}
        curated_variant_map = {v.name: v for v in all_curated_variants}
        
        for name in selected_curated_variants:
            if name not in existing_variant_names and name in curated_variant_map:
                variant_to_add = curated_variant_map[name]
                signature_variant = Variant.from_signature_variant(variant_to_add)
                
                # Register in the unified registry first
                from state import VariantSource
                AbundanceEstimatorState.register_variant(
                    name=signature_variant.name,
                    signature_mutations=signature_variant.signature_mutations,
                    source=VariantSource.CURATED
                )
    
    # ============== CUSTOM VARIANT CREATION ==============
    st.markdown("#### Compose Custom Variant")
    st.markdown("##### by selecting Signature Mutations from CovSpectrum")
    # Configure the component with compact functionality
    component_config = {
        'show_nucleotides_only': True,
        'slim_table': True,
        'show_distributions': False,
        'show_download': True,
        'show_plot': False,
        'title': "Custom Variant Composer",
        'show_title': False,
        'show_description': False,
        'default_variant': None,
        'default_min_abundance': 0.8,
        'default_min_coverage': 15
    }
    
    # Create a container for the component
    custom_container = st.container()
    
    # Render the variant signature component
    result = render_signature_composer(
        covSpectrum,
        component_config,
        session_prefix="custom_variant_",  
        container=custom_container
    )
    if result is not None:
        selected_mutations, _ = result
    else:
        selected_mutations = []


    # make button to add custom variant
    st.button("Add Custom Variant", key="add_custom_variant_button")
    # Check if the button was clicked
    if st.session_state.get("add_custom_variant_button", False):
        # Check if any mutations were selected
        if not selected_mutations:
            st.warning("Please select at least one mutation to create a custom variant.")
        else:
            # Show the selected mutations
            st.write("Selected Signature Mutations:")
            variant_query = st.session_state.get("custom_variant_variantQuery", "Custom Variant")
            
            # Check if the variant is already registered
            if AbundanceEstimatorState.is_variant_registered(variant_query):
                st.warning(f"Variant '{variant_query}' already exists in the list. Please choose a different name.")
            else:
                
                # Register directly in the variant registry
                from state import VariantSource
                AbundanceEstimatorState.register_variant(
                    name=variant_query,
                    signature_mutations=selected_mutations,
                    source=VariantSource.CUSTOM_COVSPECTRUM
                )
                
                # Add to the UI tracking for custom variants
                custom_selected = AbundanceEstimatorState.get_selected_custom_names()
                if variant_query not in custom_selected:
                    custom_selected.append(variant_query)
                    AbundanceEstimatorState.set_selected_custom_names(custom_selected)
                
                logging.info(f"Added custom variant '{variant_query}' with {len(selected_mutations)} mutations.")
                
                # Show confirmation
                st.success(f"Added custom variant '{variant_query}' with {len(selected_mutations)} mutations")
                
                # Trigger a rerun to immediately update the UI with the new variant
                st.rerun()
    
    # --- Manual Input Section ---
    with st.expander("##### Or Manual Input", expanded=False):
        manual_variant_name = st.text_input(
            "Variant Name", 
            value="", 
            max_chars=20,  # Increased max_chars slightly
            placeholder="Enter unique variant name", 
            key="manual_variant_name_input"
        )
        manual_mutations_input_str = st.text_area(
            "Signature Mutations", 
            value="", 
            placeholder="e.g., C345T, 456-, 748G", 
            help="Enter mutations separated by commas. Format: [REF]Position[ALT]. REF is optional (e.g., for insertions like 748G or deletions like 456-).",
            key="manual_mutations_input"
        )

        if st.button("Add Manual Variant", key="add_manual_variant_button"):
            # Validate variant name
            if not manual_variant_name.strip():
                st.error("Manual Variant Name cannot be empty.")
            else:
                # Process mutations
                mutations_str_list = [m.strip() for m in manual_mutations_input_str.split(',') if m.strip()]
                
                validated_signature_mutations = []
                all_mutations_valid = True

                if not mutations_str_list:
                    st.warning(f"No mutations entered for '{manual_variant_name}'. It will be added with an empty signature if the name is unique.")
                
                for mut_str in mutations_str_list:
                    # Use the improved validation method from the Mutation class
                    is_valid, error_message, mutation_data = Mutation.validate_mutation_string(mut_str)
                    
                    if is_valid:
                        validated_signature_mutations.append(mut_str) # Store original valid string
                    else:
                        st.error(f"Invalid mutation: {error_message}")
                        all_mutations_valid = False

                if all_mutations_valid:
                    # Check if the variant is already registered
                    if AbundanceEstimatorState.is_variant_registered(manual_variant_name):
                        st.warning(f"Variant '{manual_variant_name}' already exists in the list. Please choose a different name.")
                    else:
                        
                        # Register directly in the variant registry
                        from state import VariantSource
                        AbundanceEstimatorState.register_variant(
                            name=manual_variant_name,
                            signature_mutations=validated_signature_mutations,
                            source=VariantSource.CUSTOM_MANUAL
                        )
                        
                        # Add to the UI tracking for custom variants
                        custom_selected = AbundanceEstimatorState.get_selected_custom_names()
                        if manual_variant_name not in custom_selected:
                            custom_selected.append(manual_variant_name)
                            AbundanceEstimatorState.set_selected_custom_names(custom_selected)
                        
                        st.success(f"Added manual variant '{manual_variant_name}' with {len(validated_signature_mutations)} mutations.")
                        # Set flag to clear inputs on next rerun
                        AbundanceEstimatorState.clear_manual_inputs()
                        st.rerun() # Trigger rerun
    
    # ============== VARIANT VALIDATION AND PROCESSING ==============
    
    # Only show warning if combined_variants.variants is empty AND we've already loaded data
    if not combined_variants.variants and AbundanceEstimatorState.get_registered_variants():
        st.warning("Please select at least one variant from either the curated list or create a custom variant")
    
    st.markdown("---")
    
    # ============== SELECTED VARIANTS MANAGEMENT ==============
    st.subheader("Selected Variants")
    
    # URL sharing limitation info
    st.info("💡 **URL Sharing Note:** Only curated variants are included in shareable URLs. Custom variants (created via CovSpectrum or manual input) are not saved to URLs due to potential length limitations. To share complete analysis including custom variants, use the export functionality below.")
    
    # Get the current variants for display
    current_variant_names = [variant.name for variant in combined_variants.variants]
    
    if current_variant_names:
        # Create a multiselect showing current variants (user can deselect to remove)
        displayed_variant_names = st.multiselect(
            "Currently Selected Variants (Deselect to remove)",
            options=current_variant_names,
            default=current_variant_names,  # Always default to show all current variants
            help="Deselect variants to remove them from the list."
        )
        
        # Check if any variants were deselected
        variants_removed = False
        all_removed = False
        
        # Get reference to state for tracking
        curated_selected = AbundanceEstimatorState.get_selected_curated_names()
        custom_selected = AbundanceEstimatorState.get_selected_custom_names()
        
        for variant in combined_variants.variants:
            if variant.name not in displayed_variant_names:
                # Remove from unified variant registry
                AbundanceEstimatorState.unregister_variant(variant.name)
                
                # Also remove from UI tracking lists
                if variant.name in curated_selected:
                    curated_selected.remove(variant.name)
                    AbundanceEstimatorState.set_selected_curated_names(curated_selected)
                
                if variant.name in custom_selected:
                    custom_selected.remove(variant.name)
                    AbundanceEstimatorState.set_selected_custom_names(custom_selected)
                
                st.success(f"Removed variant '{variant.name}' from the list.")
                variants_removed = True
        
        # If all variants were deselected, handle this special case
        if not displayed_variant_names and current_variant_names:
            all_removed = True
            # Clear all registries
            for name in current_variant_names:
                AbundanceEstimatorState.unregister_variant(name)
            
            # Clear UI tracking lists
            AbundanceEstimatorState.set_selected_curated_names([])
            AbundanceEstimatorState.set_selected_custom_names([])
            st.warning("All variants were removed. Please select new variants from the curated list or create custom variants.")
        
        # If variants were removed, rerun to update the UI
        if variants_removed or all_removed:
            # Update URL to reflect the current curated variants selection
            current_curated = AbundanceEstimatorState.get_selected_curated_names()
            url_state.save_to_url(selected_variants=current_curated)
            st.rerun()
    else:
        st.info("No variants are currently selected. Select variants from the Curated Variant List or create Custom Variants above.")
    
    # ============== VARIANT DEBUG INFO (EXPANDABLE) ==============
    with st.expander("🔍 Variant Selection Information", expanded=False):
        st.write("**Currently Selected Variants:**")
        
        registered_variants = AbundanceEstimatorState.get_registered_variants()
        if registered_variants:
            for variant_name, variant_data in registered_variants.items():
                try:
                    col1, col2, col3 = st.columns([2, 1, 3])
                    with col1:
                        st.write(f"**{variant_name}**")
                    with col2:
                        # Display the source as a nicely formatted string
                        source_display = variant_data['source'].value.replace('_', ' ').title()
                        st.write(f"*{source_display}*")
                    with col3:
                        st.write(f"{len(variant_data['signature_mutations'])} mutations")
                except (ValueError, TypeError):
                    # Fallback for test environments where st.columns might not work properly
                    source_display = variant_data['source'].value.replace('_', ' ').title()
                    st.write(f"**{variant_name}** - *{source_display}* - {len(variant_data['signature_mutations'])} mutations")
        else:
            st.write("No custom variants registered")
        
    # ============== MUTATION-VARIANT MATRIX ==============
    st.markdown("---")
    
    if not combined_variants.variants:
        st.info("Select at least one variant to see the mutation-variant matrix and visualizations.")
        
    # Build the mutation-variant matrix
    elif combined_variants.variants:
        
        # Create the mutation-variant matrix using the utility function
        matrix_df = create_mutation_variant_matrix(combined_variants)
        
        # Create a section with two visualizations side by side
        
        # ============== VARIANT SIGNATURE COMPARISON ==============
        st.subheader("Variant Signature Comparison")
        st.write("Compare the selected variant signatures via variant matrix and venn diagrams")

        # Visualize the data in different ways
        if len(combined_variants.variants) > 1:
            
            # Create a matrix to show shared mutations between variants (triangular)
            variant_names = [variant.name for variant in combined_variants.variants]
            variant_comparison = pd.DataFrame(index=variant_names, columns=variant_names)

            # For each pair of variants, count the number of shared mutations
            # Only compute upper triangle + diagonal to avoid redundancy
            for i, variant1 in enumerate(combined_variants.variants):
                for j, variant2 in enumerate(combined_variants.variants):
                    if j >= i:  # Upper triangle + diagonal only
                        # Get the sets of mutations for each variant
                        mutations1 = set(variant1.signature_mutations)
                        mutations2 = set(variant2.signature_mutations)
                        
                        # Count number of shared mutations
                        shared_count = len(mutations1.intersection(mutations2))
                        
                        # Store in the dataframe
                        variant_comparison.iloc[i, j] = shared_count
                    else:
                        # Lower triangle: set to NaN to hide in visualization
                        variant_comparison.iloc[i, j] = np.nan
        
            # Convert numeric data to int for display, keeping NaN for lower triangle
            # Only convert the upper triangle values to int
            for i in range(len(variant_comparison.index)):
                for j in range(len(variant_comparison.columns)):
                    if j >= i and not pd.isna(variant_comparison.iloc[i, j]):
                        variant_comparison.iloc[i, j] = int(variant_comparison.iloc[i, j])
            variant_comparison_melted = variant_comparison.reset_index().melt(
                id_vars="index", 
                var_name="variant2", 
                value_name="shared_mutations"
            )
            variant_comparison_melted.columns = ["variant1", "variant2", "shared_mutations"]
            
            
            # Create two columns for the visualizations
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("#### Shared Mutations Heatmap")
                # Calculate the shared mutations for hover text (only upper triangle + diagonal)
                shared_mutations_hover = {}
                for i, variant1 in enumerate(combined_variants.variants):
                    for j, variant2 in enumerate(combined_variants.variants):
                        if j >= i:  # Only upper triangle + diagonal
                            mutations1 = set(variant1.signature_mutations)
                            mutations2 = set(variant2.signature_mutations)
                            shared = mutations1.intersection(mutations2)
                            shared_mutations_hover[(variant1.name, variant2.name)] = shared

                # Create hover text with shared mutations
                hover_text = []
                for i, variant1 in enumerate([v.name for v in combined_variants.variants]):
                    hover_row = []
                    for j, variant2 in enumerate([v.name for v in combined_variants.variants]):
                        if j >= i:  # Upper triangle + diagonal
                            count = variant_comparison.iloc[i, j]
                            shared = shared_mutations_hover.get((variant1, variant2), set())
                            
                            if variant1 == variant2:
                                text = f"<b>{variant1}</b><br>{int(count)} signature mutations"
                            else:
                                text = f"<b>{variant1} ∩ {variant2}</b><br>{int(count)} shared mutations"
                                if shared:
                                    mutations_list = list(shared)
                                    if len(mutations_list) > 10:
                                        text += f"<br>First 10 shared mutations:<br>• " + "<br>• ".join(mutations_list[:10]) + f"<br>...and {len(mutations_list)-10} more"
                                    else:
                                        text += "<br>Shared mutations:<br>• " + "<br>• ".join(mutations_list)
                        else:
                            # Lower triangle: empty hover text
                            text = ""
                        
                        hover_row.append(text)
                    hover_text.append(hover_row)

                # Get min and max values for better color mapping (excluding NaN)
                valid_values = variant_comparison.values[~pd.isna(variant_comparison.values)]
                min_val = valid_values.min()
                max_val = valid_values.max()
                
                # Create annotation text with adaptive text color (only for upper triangle + diagonal)
                annotations = []
                for i in range(len(variant_comparison.index)):
                    for j in range(len(variant_comparison.columns)):
                        if j >= i and not pd.isna(variant_comparison.iloc[i, j]):  # Upper triangle + diagonal
                            value = variant_comparison.iloc[i, j]
                            # Normalize value between 0 and 1
                            normalized_val = (value - min_val) / (max_val - min_val) if max_val > min_val else 0
                            # Adjust threshold based on the Blues colorscale - text is white above this normalized value
                            text_color = "white" if normalized_val > 0.5 else "black"
                            
                            annotations.append(
                                dict(
                                    x=j,
                                    y=i,
                                    text=str(int(value)),
                                    showarrow=False,
                                    font=dict(color=text_color, size=14)
                                )
                            )

                # Determine size based on number of variants (square plot)
                size = max(350, min(500, 100 * len(combined_variants.variants)))

                # Create Plotly heatmap (triangular)
                fig = go.Figure(data=go.Heatmap(
                    z=variant_comparison.values,
                    x=variant_comparison.columns,
                    y=variant_comparison.index,
                    colorscale='Blues',
                    text=hover_text,
                    hoverinfo='text',
                    showscale=False,  # Remove colorbar/legend
                    hoverongaps=False  # Don't show hover for NaN values
                ))

                # Update layout
                fig.update_layout(
                    xaxis=dict(title='Variant', showgrid=False),
                    yaxis=dict(title='Variant', showgrid=False),
                    height=size,
                    width=size,
                    margin=dict(l=50, r=20, t=50, b=20),
                    annotations=annotations
                )

                # Display the interactive Plotly chart in Streamlit
                st.plotly_chart(fig)
            
            # Venn Diagram in the second column (supports 2-3 variants)
            with col2:
                if 2 <= len(combined_variants.variants) <= 3:
                    st.markdown("#### Mutation Overlap")
                    
                    # Matplotlib is already imported at the top
                    matplotlib.use('agg')  # Set non-interactive backend
                    
                    # Set a professional style for the plots
                    plt.style.use('seaborn-v0_8-whitegrid')  # Modern, clean style
                    
                    if len(combined_variants.variants) == 2:
                        from matplotlib_venn import venn2
                        
                        # Create sets of mutations for each variant
                        set1 = set(combined_variants.variants[0].signature_mutations)
                        set2 = set(combined_variants.variants[1].signature_mutations)
            
                        # Very small fixed size - try even smaller
                        fig_venn, ax_venn = plt.subplots(figsize=(3, 3), dpi=100)
                        
                        # Fix the typing issue by explicitly creating a tuple of exactly 2 elements
                        variant_name1 = combined_variants.variants[0].name
                        variant_name2 = combined_variants.variants[1].name
                        variant_labels = (variant_name1, variant_name2)
                        
                        venn2((set1, set2), variant_labels, ax=ax_venn)
                        
                        # Adjust layout to be more compact
                        plt.tight_layout(pad=1.5)
                        
                        # Add a light gray border
                        for spine in ax_venn.spines.values():
                            spine.set_visible(True)
                            spine.set_color('#f0f0f0')
                        
                        # Display the venn diagram with controlled size using columns
                        venn_col1, venn_col2, venn_col3 = st.columns([1, 2, 1])
                        with venn_col2:
                            st.pyplot(fig_venn)
                        
                    elif len(combined_variants.variants) == 3:
                        from matplotlib_venn import venn3
                        
                        # Create sets of mutations for each variant - extract exactly 3 sets as required by venn3
                        set1 = set(combined_variants.variants[0].signature_mutations)
                        set2 = set(combined_variants.variants[1].signature_mutations)
                        set3 = set(combined_variants.variants[2].signature_mutations)
                        
                        # Very small fixed size - try even smaller
                        fig_venn, ax_venn = plt.subplots(figsize=(3, 3), dpi=100)
                        
                        # Fix the typing issue by explicitly creating a tuple of exactly 3 elements
                        variant_name1 = combined_variants.variants[0].name
                        variant_name2 = combined_variants.variants[1].name
                        variant_name3 = combined_variants.variants[2].name
                        variant_labels = (variant_name1, variant_name2, variant_name3)
                        
                        venn3((set1, set2, set3), variant_labels, ax=ax_venn)
                        
                        # Adjust layout to be more compact
                        plt.tight_layout(pad=1.5)
                        
                        # Add a light gray border
                        for spine in ax_venn.spines.values():
                            spine.set_visible(True)
                            spine.set_color('#f0f0f0')
                        
                        # Display the venn diagram with controlled size using columns
                        venn_col1, venn_col2, venn_col3 = st.columns([1, 2, 1])
                        with venn_col2:
                            st.pyplot(fig_venn)
                else:
                    st.markdown("#### Mutation Overlap")
                    st.info("Venn diagram is only available for 2-3 variants")
            

            # 3. Mutation-Variant Matrix Visualization (heatmap) - Collapsible
            with st.expander("Variant-Signatures Bitmap Visualization", expanded=False):

                st.write("This heatmap shows which mutations (rows) are present in each variant (columns). Blue cells indicate the mutation is present.")
                
                # First prepare the data in a suitable format
                binary_matrix = matrix_df.set_index("Mutation")
        
                # Use Plotly for a more interactive visualization
                fig = go.Figure(data=go.Heatmap(
                    z=binary_matrix.values,
                    x=binary_matrix.columns,
                    y=binary_matrix.index,
                    colorscale=[[0, 'white'], [1, '#1E88E5']],  # Match the color scheme
                    showscale=False,  # Hide color scale bar
                    hoverongaps=False
                ))

                # Calculate dimensions based on data size
                num_mutations = len(binary_matrix.index)
                num_variants = len(binary_matrix.columns)

                # Customize layout with settings to ensure all labels are visible
                fig.update_layout(
                    title='Mutation-Variant Matrix',
                    xaxis=dict(
                        title='Variant',
                        side='top',  # Show x-axis on top
                    ),
                    yaxis=dict(
                        title='Mutation',
                        automargin=True,  # Automatically adjust margins for labels
                        tickmode='array',  # Force all ticks
                        tickvals=list(range(len(binary_matrix.index))),  # Positions for each mutation
                        ticktext=binary_matrix.index,  # Actual mutation labels
                    ),
                    height=max(500, min(2000, 25 * num_mutations)),  # Dynamic height based on mutations
                    width=max(600, 120 * num_variants),  # Dynamic width based on variants
                    margin=dict(l=150, r=20, t=50, b=20),  # Increase left margin for y labels
                )
                
                # Add custom hover text
                hover_text = []
                for i, mutation in enumerate(binary_matrix.index):
                    row_hover = []
                    for j, variant in enumerate(binary_matrix.columns):
                        if binary_matrix.iloc[i, j] == 1:
                            text = f"Mutation: {mutation}<br>Variant: {variant}<br>Status: Present"
                        else:
                            text = f"Mutation: {mutation}<br>Variant: {variant}<br>Status: Absent"
                        row_hover.append(text)
                    hover_text.append(row_hover)
                
                fig.update_traces(hoverinfo='text', text=hover_text)
                
                # Display the interactive Plotly chart in Streamlit
                st.plotly_chart(fig)
        else:
            st.warning("At least two variants are required to visualize the mutation-variant matrix.")
    
        # ============== EXPORT FUNCTIONALITY ==============
        # Export functionality
        st.subheader("Export Variant Signatures")
        
        if combined_variants.variants:
            # Convert to CSV for download
            csv = matrix_df.to_csv(index=False)
            st.download_button(
                label="Download Mutation-Variant Matrix (CSV)",
                data=csv,
                file_name="mutation_variant_matrix.csv",
                mime="text/csv",
            )
        else:
            st.info("Select at least one variant to enable export functionality.")
        
        st.markdown("---")

        # ============== ANALYSIS CONFIGURATION ==============
        st.subheader("Analysis Configuration & Execution")
        st.write("Configure the parameters for data fetching and variant abundance estimation, then run the complete analysis.")
        
        # Add information about the process
        st.info("💡 This analysis will: 1) Fetch mutation counts and coverage data from Loculus, 2) Estimate variant abundances using LolliPop deconvolution")
        
        # Initialize the API client first
        server_ip = get_wiseloculus_url()
        wiseLoculus = WiseLoculusLapis(server_ip)
        
        # Configuration sections in columns
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("#### 📊 Data Parameters")
            
            # Date range input
            # Get dynamic date range from API with bounds to enforce limits
            default_start, default_end, min_date, max_date = wiseLoculus.get_cached_date_range_with_bounds("abundance")
            
            date_range = st.date_input(
                "Date Range",
                value=[default_start, default_end],
                min_value=min_date,
                max_value=max_date,
                help="Select the date range for fetching mutation counts and coverage data from Loculus."
            )

            locations = wiseLoculus.fetch_locations()
            if locations is None:
                st.error("Unable to fetch locations. Please check your connection to the Loculus server.")
                selected_locations = []
            else:
                selected_locations = st.multiselect(
                    "Locations", 
                    options=locations,
                    default=locations,  # Default to all locations selected
                    help="Select one or more sampling locations for data analysis. Multiple locations will be processed together."
                )
                
                # Show warning if more than one location is selected
                if len(selected_locations) > 1:
                    st.warning(f"⚠️ Multiple locations selected ({len(selected_locations)} locations). Multi-location analysis is currently experimental and may take longer to process.")
                elif len(selected_locations) == 0:
                    st.error("Please select at least one location for analysis.")
                else:
                    st.success(f"✅ Single location selected: {selected_locations[0]}")
            
            # For backward compatibility, set location to the first selected location or None
            location = selected_locations[0] if selected_locations else None
        
        with col2:
            st.markdown("#### ⚙️ Deconvolution Parameters")
            
            # Bootstrap iterations with preset options
            bootstrap_options = {
                "Rapid (Fast)": 50,
                "Standard": 100, 
                "Reliable (Slower)": 300
            }
            
            selected_option = st.radio(
                "Bootstrap Iterations",
                options=list(bootstrap_options.keys()),
                index=1,  # Default to "Standard" (100 bootstraps)
                help="Choose the number of bootstrap iterations for confidence intervals. More iterations provide better estimates but take longer to compute."
            )
            
            bootstraps = bootstrap_options[selected_option]

            # Bandwidth parameter for gaussian kernel smoothing
            bandwidth_options = {
                "Narrow": 10,
                "Medium": 20,  
                "Wide": 30,
            }
            
            selected_bandwidth_option = st.radio(
                "Bandwidth (Gaussian Kernel Smoothing)",
                options=list(bandwidth_options.keys()),
                index=0,  # Default to "Narrow" (10)
                help="Controls the smoothing applied to time series data. The optimal smoothing is highly dependent on the noise in the data. Low noise allows for narrower timeranges, while higher bandwidths help counter noise. Information is shared across time points. Narrow bandwidth preserves short-term variations (1-2 months), wide bandwidth smooths long-term trends (3+ months). Choose based on your timeframe and number of data points."
            )
            
            bandwidth = bandwidth_options[selected_bandwidth_option]
            
            # Show the actual bandwidth value selected
            st.caption(f"Selected: Bandwidth = {bandwidth}")
            
            # Additional guidance based on date range
            if len(date_range) == 2:
                days_diff = (date_range[1] - date_range[0]).days
                if days_diff <= 60 and bandwidth > 10:
                    st.info("💡 For timeframes ≤2 months, consider using 'Narrow' bandwidth for better short-term variation capture.")
                elif days_diff > 90 and bandwidth < 20:
                    st.info("💡 For timeframes >3 months, consider using 'Wide' bandwidth for better trend smoothing.")


            
            # Show the actual number selected
            st.caption(f"Selected: {bootstraps} bootstrap iterations")
            
            st.markdown("**Method:** LolliPop Deconvolution")
            st.caption("Kernel-based deconvolution that leverages time series data to generate high-confidence abundance estimates despite noise in wastewater samples.")
        
        # Multi-location session state is already initialized above

        # Initialize multi-location session state variables
        if 'location_data' not in st.session_state:
            st.session_state.location_data = {}
        
        if 'location_tasks' not in st.session_state:
            st.session_state.location_tasks = {}  # {location: task_id}
        
        if 'location_results' not in st.session_state:
            st.session_state.location_results = {}  # {location: result_data}

        # Single action button
        st.markdown("---")
        
        # Button logic depends on whether we already have results
        if st.session_state.location_results:
            # If we have results, show a "Start New Analysis" button instead
            if st.button("🔄 Start New Analysis", help="Clear current results and start a new complete analysis", type="primary"):
                # Clear all analysis state to reset the workflow
                st.session_state.location_data = {}
                st.session_state.location_tasks = {}
                st.session_state.location_results = {}
                
                # Also clear any data hash to force recomputation
                if 'last_data_hash' in st.session_state:
                    del st.session_state.last_data_hash
                st.rerun()  # Rerun to show the parameters and Run Analysis button
        else:
            # Only show Run Analysis button if we don't have results yet
            if st.button("Run Complete Analysis", help="Fetch data and estimate variant abundances", type="primary"):
                # Validate inputs
                if not selected_locations:
                    st.error("Please select at least one location. Unable to fetch data without selecting a location.")
                    st.stop()
                
                if len(date_range) != 2:
                    st.error("Please select a valid date range with both start and end dates.")
                    st.stop()
                
                # Show information about the analysis scope
                if len(selected_locations) > 1:
                    st.info(f"Running analysis for {len(selected_locations)} locations: {', '.join(selected_locations)}")
                else:
                    st.info(f"Running analysis for location: {selected_locations[0]}")
                
                # Get the latest mutation list
                mutations = matrix_df["Mutation"].tolist()
                
                # Convert date range to datetime tuples
                start_date = datetime.combine(date_range[0], datetime.min.time())
                end_date = datetime.combine(date_range[1], datetime.min.time())
                datetime_range = (start_date, end_date)
                
                # Step 1: Fetch data
                with st.spinner('Step 1/2: Fetching mutation counts and coverage data from Loculus...'):
                    try:
                        # Import the new multi-location utility
                        from utils.multi_location import fetch_multi_location_data
                        
                        # Fetch data for all locations (parallel or single)
                        st.session_state.location_data = asyncio.run(fetch_multi_location_data(
                            wiseLoculus,
                            mutations,
                            MutationType.NUCLEOTIDE,  
                            datetime_range,
                            selected_locations,
                            interval="daily"
                        ))
                        
                        # Data successfully fetched for all locations
                        
                        st.success("✅ Data fetched successfully for all locations!")
                            
                    except Exception as e:
                        st.error(f"❌ Error fetching data: {str(e)}")
                        st.stop()
                
                # Step 2: Run deconvolution for each location
                if st.session_state.location_data:
                    with st.spinner('Step 2/2: Starting variant abundance estimation for all locations...'):
                        try:
                            # Submit separate tasks for each location
                            st.session_state.location_tasks = {}
                            
                            for location, location_counts_df in st.session_state.location_data.items():
                                st.write(f"Starting deconvolution for {location}...")
                                
                                # Convert DataFrames to base64-encoded pickle strings for serialization
                                counts_pickle = base64.b64encode(pickle.dumps(location_counts_df)).decode('utf-8')
                                matrix_pickle = base64.b64encode(pickle.dumps(matrix_df)).decode('utf-8')
                                
                                # Submit the task to Celery worker
                                task = celery_app.send_task(
                                    'tasks.run_deconvolve',
                                    kwargs={
                                        'mutation_counts_df': counts_pickle,
                                        'mutation_variant_matrix_df': matrix_pickle,
                                        'bootstraps': bootstraps,
                                        'bandwidth': bandwidth,  # Add bandwidth parameter
                                        'locationName': location  # Add location name
                                    }
                                )
                                
                                # Store task ID for this location
                                st.session_state.location_tasks[location] = task.id
                                st.success(f"✅ Analysis started for {location}!")
                            
                            # Task submission completed for all locations
                            
                        except Exception as e:
                            st.error(f"❌ Error starting deconvolution: {str(e)}")
                            st.stop()
        
        # Show download options for fetched data (if available)
        if st.session_state.location_data:
            with st.expander("📥 Download Raw Data", expanded=False):
                st.write("Download the fetched mutation counts and coverage data:")
                
                if len(st.session_state.location_data) == 1:
                    # Single location - show simple download
                    locationName = list(st.session_state.location_data.keys())[0]
                    location_data = st.session_state.location_data[locationName]
                    
                    # Create columns for download buttons
                    col1, col2 = st.columns(2)
                    
                    # 1. CSV Download
                    with col1:
                        # Make sure to preserve index for dates and mutations
                        csv = location_data.to_csv(index=True)
                        st.download_button(
                            label=f"Download {locationName} as CSV",
                            data=csv,
                            file_name=f'mutation_counts_coverage_{locationName.replace(" ", "_")}.csv',
                            mime='text/csv',
                            help="Download all data as a single CSV file with preserved indices."
                        )
                    
                    # 2. JSON Download
                    with col2:
                        # Convert to JSON structure - using 'split' format to preserve indices
                        json_data = location_data.to_json(orient='split', date_format='iso', index=True)
                        
                        st.download_button(
                            label=f"Download {locationName} as JSON",
                            data=json_data,
                            file_name=f'mutation_counts_coverage_{locationName.replace(" ", "_")}.json',
                            mime='application/json',
                            help="Download data as a JSON file that preserves dates and mutation indices."
                        )
                else:
                    # Multiple locations - show combined and individual downloads
                    st.markdown("#### Combined Download (All Locations)")
                    
                    # Combine all location data
                    combined_dfs = []
                    total_data_points = 0
                    
                    for locationName, location_data in st.session_state.location_data.items():
                        try:
                            if location_data is not None and not location_data.empty:
                                # Add location column to each dataframe
                                location_df = location_data.copy()
                                location_df['location'] = locationName
                                combined_dfs.append(location_df)
                                total_data_points += len(location_df)
                                st.caption(f"✅ {locationName}: {len(location_df)} data points")
                            else:
                                st.caption(f"⚠️ {locationName}: No data available")
                        except Exception as e:
                            st.caption(f"❌ {locationName}: Error processing data - {str(e)}")
                    
                    if combined_dfs:
                        # Concatenate all dataframes
                        try:
                            combined_data = pd.concat(combined_dfs, ignore_index=False)
                            st.success(f"📊 Combined {len(combined_dfs)} locations with {total_data_points} total data points")
                            
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                combined_csv = combined_data.to_csv(index=True)
                                st.download_button(
                                    label="📄 Download All Locations (CSV)",
                                    data=combined_csv,
                                    file_name='mutation_counts_coverage_all_locations.csv',
                                    mime='text/csv',
                                    help="Download combined data from all locations with location column"
                                )
                            
                            with col2:
                                combined_json = combined_data.to_json(orient='split', date_format='iso', index=True)
                                st.download_button(
                                    label="📋 Download All Locations (JSON)",
                                    data=combined_json,
                                    file_name='mutation_counts_coverage_all_locations.json',
                                    mime='application/json',
                                    help="Download combined data from all locations as JSON"
                                )
                        except Exception as e:
                            st.error(f"Error combining data: {str(e)}")
                            st.warning("Unable to create combined download. Check debug information above.")
                    else:
                        st.warning("No valid data found for combined download. Check individual locations in debug section above.")
        
        st.markdown("---")
        
        # ============== RESULTS SECTION ==============
        st.subheader("Analysis Results")

        # Check if we have any location tasks
        if st.session_state.location_tasks:
            # Check processing status first and display it prominently
            incomplete_tasks = []
            for location, task_id in st.session_state.location_tasks.items():
                if location not in st.session_state.location_results:
                    try:
                        task = celery_app.AsyncResult(task_id)
                        if not task.ready():
                            incomplete_tasks.append(location)
                    except Exception:
                        incomplete_tasks.append(location)  # Consider as incomplete if we can't check
            
            if incomplete_tasks:
                st.info(f"⏳ Still processing: {', '.join(incomplete_tasks)}")
                
                # Auto-refresh every 5 seconds if there are incomplete tasks
                from streamlit_autorefresh import st_autorefresh
                st_autorefresh(interval=5000, key="multi_location_autorefresh")
            else:
                st.success("🎉 All locations have completed analysis!")
                        
            # Import and render the multi-location results component
            from components.multi_location_results import render_location_results_tabs
            
            render_location_results_tabs(
                st.session_state.location_tasks,
                st.session_state.location_results,
                celery_app,
                redis_client
            )

            # ============== COVVFIT FITNESS INFERENCE SECTION ==============
            # covvfit needs completed deconvolution for all locations.
            all_locations_done = (
                    len(st.session_state.location_results) == len(st.session_state.location_tasks)
                    and len(st.session_state.location_results) > 0
            )

            if all_locations_done:
                st.markdown("---")
                st.subheader("Variant fitness inference (covvfit)")

                # Initialize covvfit session state
                if 'covvfit_task_id' not in st.session_state:
                    st.session_state.covvfit_task_id = None
                if 'covvfit_result' not in st.session_state:
                    st.session_state.covvfit_result = None

                col_a, col_b = st.columns(2)
                with col_a:
                    max_days = st.slider(
                        "Days of past data to use",
                        min_value=60, max_value=365, value=180, step=10,
                        help="How many days of historical data covvfit fits the model to.",
                    )
                with col_b:
                    horizon = st.slider(
                        "Days to forecast",
                        min_value=30, max_value=180, value=90, step=10,
                        help="How many days into the future covvfit projects.",
                    )

                if st.button("Run fitness inference", type="primary"):
                    matrix_pickle = base64.b64encode(pickle.dumps(matrix_df)).decode("utf-8")

                    location_data_pickled = {}
                    for loc, counts_df in st.session_state.location_data.items():
                        location_data_pickled[loc] = base64.b64encode(
                            pickle.dumps(counts_df)
                        ).decode("utf-8")

                    task = celery_app.send_task(
                        "tasks.run_covvfit",
                        kwargs={
                            "location_data": location_data_pickled,
                            "matrix_df": matrix_pickle,
                            "max_days": max_days,
                            "horizon": horizon,
                        },
                    )

                    st.session_state.covvfit_task_id = task.id
                    st.session_state.covvfit_result = None
                    st.rerun()

                # Poll / display covvfit result
                if st.session_state.covvfit_task_id is not None:
                    if st.session_state.covvfit_result is not None:
                        # We already have the result — display it
                        # TODO: with a single location covvfit places the legend far
                        # from the panel, leaving a large gap. Could be tuned via a
                        # covvfit plot config (right / panel_width / wspace). Left as-is
                        # for now since multi-location runs (the common case) look fine.
                        figure_png_b64 = st.session_state.covvfit_result["figure_png"]
                        st.image(base64.b64decode(figure_png_b64))

                        st.markdown("#### Download fitness data")

                        if "figure_pdf" in st.session_state.covvfit_result:
                            st.download_button(
                                label="Download figure (PDF)",
                                data=base64.b64decode(st.session_state.covvfit_result["figure_pdf"]),
                                file_name="covvfit_figure.pdf",
                                mime="application/pdf",
                            )

                        if "pairwise_fitnesses_csv" in st.session_state.covvfit_result:
                            st.download_button(
                                label="Download fitness results (CSV)",
                                data=st.session_state.covvfit_result["pairwise_fitnesses_csv"],
                                file_name="covvfit_pairwise_fitnesses.csv",
                                mime="text/csv",
                            )

                    else:
                        # Task is running — check status
                        covvfit_task = celery_app.AsyncResult(st.session_state.covvfit_task_id)
                        if covvfit_task.ready():
                            try:
                                st.session_state.covvfit_result = covvfit_task.get()
                                st.rerun()
                            except Exception as e:
                                st.error(f"covvfit failed: {str(e)}")
                        else:
                            # Show progress from Redis
                            progress_raw = redis_client.get(
                                f"task_progress:{st.session_state.covvfit_task_id}"
                            )
                            if progress_raw:
                                progress = json.loads(progress_raw)
                                st.info(f"⏳ {progress.get('status', 'Running covvfit...')}")
                            else:
                                st.info("⏳ Running covvfit...")

                            from streamlit_autorefresh import st_autorefresh
                            st_autorefresh(interval=5000, key="covvfit_autorefresh")

        else:
            # No tasks have been started yet
            if not st.session_state.location_data:
                st.info("Run the complete analysis above to see results here.")
            else:
                st.info("Data has been fetched. Run the complete analysis to estimate variant abundances.")
            
if __name__ == "__main__":
    app()