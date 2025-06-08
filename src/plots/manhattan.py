"""
Manhattan plot generation functions for STRXplorer
"""
import pandas as pd
import plotly.graph_objects as go
import numpy as np


def create_manhattan_plot(
    data: pd.DataFrame,
    trait_name: str,
    target_chrom: str,
    target_pos: int,
    window_size: int,
):
    """Create Manhattan plot from database data"""

    fig = go.Figure()

    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size

    if not data.empty:
        # Create hover text
        hover_text = []
        for _, row in data.iterrows():
            text = f"Chr{row['chrom']}:{row['pos']}<br>"
            text += f"P-value: {row['p_value']:.2e}<br>"
            text += f"-log10(p): {row['neg_log_p']:.2f}<br>"
            text += f"ID: {row['variant_id']}"
            if pd.notna(row["beta"]):
                text += f"<br>Beta: {row['beta']:.3f}"
            hover_text.append(text)

        # Add scatter plot
        fig.add_trace(
            go.Scatter(
                x=data["pos"],
                y=data["neg_log_p"],
                mode="markers",
                marker=dict(
                    color="rgba(31, 119, 180, 0.7)", size=6, line=dict(width=0)
                ),
                name="Variants",
                text=hover_text,
                hoverinfo="text",
                hovertemplate="%{text}<extra></extra>",
            )
        )

        # Add vertical line at target position
        max_y = data["neg_log_p"].max() * 1.05 if len(data) > 0 else 10
        fig.add_shape(
            type="line",
            x0=target_pos,
            x1=target_pos,
            y0=0,
            y1=max_y,
            line=dict(color="red", width=2, dash="dash"),
        )

        # Add significance threshold
        significance_line = -np.log10(5e-8)
        fig.add_shape(
            type="line",
            x0=start_pos,
            x1=end_pos,
            y0=significance_line,
            y1=significance_line,
            line=dict(color="blue", width=1.5, dash="dash"),
        )

        # Add annotation for significance threshold
        fig.add_annotation(
            x=start_pos + (window_size * 0.05),
            y=significance_line,
            text="p = 5e-8",
            showarrow=False,
            yshift=10,
            font=dict(color="blue", size=12),
        )

    # Update layout
    fig.update_layout(
        title=f'Manhattan Plot: {trait_name.replace("_", " ").title()} ±{window_size/1000:.0f}kb',
        xaxis=dict(
            title=f"Position on Chr{target_chrom} (bp)",
            range=[start_pos, end_pos],
            showgrid=True,
            gridcolor="lightgray",
        ),
        yaxis=dict(title="-log₁₀(p-value)", showgrid=True, gridcolor="lightgray"),
        height=500,
        hovermode="closest",
        plot_bgcolor="white",
        showlegend=False,
    )

    return fig


def create_mini_manhattan_plot(data, trait_name, chrom, target_pos, repeat_id):
    """Create a smaller Manhattan plot for the trait overview grid"""
    fig = go.Figure()

    if not data.empty:
        # Create hover text
        hover_text = []
        for _, row in data.iterrows():
            text = f"Chr{row['chrom']}:{row['pos']}<br>"
            text += f"P: {row['p_value']:.2e}<br>"
            text += f"-log10(p): {row['neg_log_p']:.2f}"
            hover_text.append(text)

        # Add scatter plot
        fig.add_trace(
            go.Scatter(
                x=data["pos"],
                y=data["neg_log_p"],
                mode="markers",
                marker=dict(
                    color="rgba(31, 119, 180, 0.6)", size=3, line=dict(width=0)
                ),
                text=hover_text,
                hoverinfo="text",
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            )
        )

        # Add vertical line at target position
        max_y = data["neg_log_p"].max() * 1.05 if len(data) > 0 else 10
        fig.add_shape(
            type="line",
            x0=target_pos,
            x1=target_pos,
            y0=0,
            y1=max_y,
            line=dict(color="red", width=1, dash="dash"),
        )

        # Add significance threshold
        significance_line = -np.log10(5e-8)
        fig.add_shape(
            type="line",
            x0=data["pos"].min(),
            x1=data["pos"].max(),
            y0=significance_line,
            y1=significance_line,
            line=dict(color="blue", width=1, dash="dash"),
        )

    # Compact layout for grid display
    fig.update_layout(
        title=f"Chr{chrom}:{target_pos:,}",
        title_font_size=12,
        xaxis=dict(
            title="Position (bp)",
            title_font_size=10,
            tickfont_size=8,
            showgrid=True,
            gridcolor="lightgray",
        ),
        yaxis=dict(
            title="-log₁₀(p)",
            title_font_size=10,
            tickfont_size=8,
            showgrid=True,
            gridcolor="lightgray",
        ),
        height=300,
        width=400,
        margin=dict(l=50, r=20, t=40, b=40),
        plot_bgcolor="white",
        hovermode="closest",
    )

    return fig
