#!/usr/bin/env python3
#plotly_rewrite.py
import argparse
import sqlite3
import json
import plotly.graph_objects as go
import numpy as np


def parse_float_or_nan(x):
    """Convert string values to float or NaN."""
    if isinstance(x, str) and x.lower() == "nan":
        return float("nan")
    return float(x)


import sqlite3
import json

def query_allele_data(db_path, repeat_id):
    """
    Queries the SQLite database for allele data using repeat_id.
    Returns:
      - dosage_dict: counts per summed length
      - mean_dict: mean phenotype values per summed length
      - ci_dict: confidence intervals per summed length
      - phenotype: the phenotype string stored in the table
      - trait_name: the trait_name string stored in the table
    """
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # pull data_json plus phenotype & trait_name
    cur.execute(
        """
        SELECT data_json, phenotype, trait_name
          FROM locus_data
         WHERE repeat_id = ?
        """,
        (repeat_id,)
    )
    row = cur.fetchone()
    conn.close()

    if not row:
        print(f"[WARNING] No data found for repeat_id: {repeat_id}")
        return None, None, None, None, None

    data_json, phenotype, trait_name = row
    data = json.loads(data_json)
    # -- your existing logic --
    dosage_dict = json.loads(data.get("sample_count_per_summed_length", "{}"))
    mean_col_name = next((c for c in data.keys() if c.startswith("mean_")), None)
    mean_dict = json.loads(data.get(mean_col_name, "{}")) if mean_col_name else {}
    ci_dict = json.loads(data.get("summed_length_0.05_alpha_CI", "{}"))

    return dosage_dict, mean_dict, ci_dict, phenotype, trait_name



def generate_figure_plotly(
    dosage_dict, mean_dict, ci_dict, trait_name, ci_style="error_bars", color_by_samples=False
):
    """
    Creates a Plotly figure with flexible visualization options.

    Parameters:
    -----------
    dosage_dict : dict
        Dictionary mapping allele lengths to sample counts
    mean_dict : dict
        Dictionary mapping allele lengths to mean phenotype values
    ci_dict : dict
        Dictionary mapping allele lengths to confidence intervals
    ci_style : str
        Confidence interval style: "error_bars" or "filled"
    color_by_samples : bool
        Whether to color markers by sample size
    """
    sorted_alleles = sorted(dosage_dict.keys(), key=float)
    if not sorted_alleles:
        print("[ERROR] No valid allele data found.")
        return None

    ci_lower = [ci_dict.get(str(a), [None, None])[0] for a in sorted_alleles]
    ci_upper = [ci_dict.get(str(a), [None, None])[1] for a in sorted_alleles]
    mean_vals = [mean_dict.get(str(a), None) for a in sorted_alleles]

    # Convert sample counts to numeric values
    sample_counts = [float(dosage_dict.get(str(a))) for a in sorted_alleles]

    fig = go.Figure()

    # Create hover text with count information
    hover_text = []
    for i, allele in enumerate(sorted_alleles):
        text = f"Length: {allele}<br>"
        text += f"Mean: {mean_vals[i]:.2f}<br>"
        text += f"95% CI: [{ci_lower[i]:.2f}, {ci_upper[i]:.2f}]<br>"
        text += f"Sample count: {sample_counts[i]:.0f}"
        hover_text.append(text)

    # Determine marker color and size based on sample size option
    marker_color = "black"
    marker_size = 8

    if color_by_samples and len(sample_counts) > 0 and max(sample_counts) > 0:
        # Normalize for visual encoding
        normalized_counts = np.array(sample_counts) / max(sample_counts)

        # Map to size and color intensity
        marker_size = [8 + nc * 12 for nc in normalized_counts]  # Size from 8 to 20

        # Use blue coloring with varying opacity
        marker_color = [
            f"rgba(31, 119, 180, {0.4 + nc * 0.6})" for nc in normalized_counts
        ]
        line_color = "rgba(31, 119, 180, 0.7)"
    else:
        # Fixed black color
        marker_color = "black"
        line_color = "black"

    # Apply different CI visualization based on style
    if ci_style == "error_bars":
        # Use error bars
        error_y = [ci_upper[i] - mean_vals[i] for i in range(len(sorted_alleles))]
        error_y_minus = [mean_vals[i] - ci_lower[i] for i in range(len(sorted_alleles))]

        error_color = "rgba(31, 119, 180, 0.3)" if color_by_samples else "black"

        fig.add_trace(
            go.Scatter(
                x=sorted_alleles,
                y=mean_vals,
                mode="lines+markers",
                error_y=dict(
                    type="data",
                    array=error_y,
                    arrayminus=error_y_minus,
                    visible=True,
                    color=error_color,
                ),
                line=dict(color=line_color, width=2),
                marker=dict(
                    color=marker_color,
                    size=marker_size,
                    line=dict(
                        width=1, color="white" if color_by_samples else line_color
                    ),
                ),
                name="95% CI",
                text=hover_text,
                hoverinfo="text",
            )
        )

    elif ci_style == "filled":
        # Fill between CI
        # Upper CI line (invisible)
        fig.add_trace(
            go.Scatter(
                x=sorted_alleles,
                y=ci_upper,
                mode="lines",
                line=dict(color="rgba(255, 0, 0, 0)"),
                showlegend=False,
                hoverinfo="none",
            )
        )

        # Lower CI line with fill between
        fig.add_trace(
            go.Scatter(
                x=sorted_alleles,
                y=ci_lower,
                mode="lines",
                line=dict(color="rgba(255, 0, 0, 0)"),
                fill="tonexty",
                fillcolor="rgba(255, 0, 0, 0.2)",
                name="95% CI",
                hoverinfo="none",
            )
        )

        # Mean line with markers
        fig.add_trace(
            go.Scatter(
                x=sorted_alleles,
                y=mean_vals,
                mode="lines+markers",
                line=dict(color=line_color, width=2),
                marker=dict(
                    color=marker_color,
                    size=marker_size,
                    line=dict(
                        width=1, color="white" if color_by_samples else line_color
                    ),
                ),
                name="Mean",
                text=hover_text,
                hoverinfo="text",
            )
        )

    # Add color scale legend if coloring by sample size
    if color_by_samples and len(sample_counts) > 0 and max(sample_counts) > 0:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    size=0,
                    colorscale=[
                        [0, "rgba(31, 119, 180, 0.4)"],
                        [1, "rgba(31, 119, 180, 1.0)"],
                    ],
                    colorbar=dict(
                        title="Sample Count",
                        thickness=15,
                        len=0.5,
                        y=0.5,
                        yanchor="middle",
                        titleside="right",
                    ),
                    cmin=0,
                    cmax=max(sample_counts),
                    showscale=True,
                ),
                hoverinfo="none",
                showlegend=False,
            )
        )

    # Enhanced layout
    fig.update_layout(
        xaxis_title="Sum of allele lengths (repeat copies)",
        yaxis_title=trait_name,
        showlegend=True,
        plot_bgcolor="white",
        xaxis=dict(
            showgrid=True,
            gridwidth=0.5,
            gridcolor="lightgray",
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor="black",
            tickformat=".2f"
        ),
        yaxis=dict(
            showgrid=True,
            gridwidth=0.5,
            gridcolor="lightgray",
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor="black",
        ),
    )

    return fig


def filter_allele_data(
    dosage_dict,
    mean_dict,
    ci_dict,
    count_threshold=100,
    max_relative_ci_range=None,
    min_length=None,
    max_length=None,
):
    """
    Enhanced filtering with length range options.
    """
    filtered_dosage = {}
    filtered_mean = {}
    filtered_ci = {}

    for allele in dosage_dict.keys():
        # Apply length range filter
        allele_val = float(allele)
        if min_length is not None and allele_val < min_length:
            continue
        if max_length is not None and allele_val > max_length:
            continue

        count_val = parse_float_or_nan(dosage_dict[allele])
        if count_val < count_threshold:
            continue

        lower = parse_float_or_nan(ci_dict.get(allele, [None, None])[0])
        upper = parse_float_or_nan(ci_dict.get(allele, [None, None])[1])
        mean_val = parse_float_or_nan(mean_dict.get(allele, None))

        if np.isnan(lower) or np.isnan(upper) or np.isnan(mean_val):
            continue

        if (
            max_relative_ci_range is not None
            and ((upper - lower) / mean_val) > max_relative_ci_range
        ):
            continue

        filtered_dosage[allele] = dosage_dict[allele]
        filtered_mean[allele] = mean_dict[allele]
        filtered_ci[allele] = ci_dict[allele]

    return filtered_dosage, filtered_mean, filtered_ci
