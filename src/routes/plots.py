"""
Plot route handlers for STRXplorer
"""
import os
from flask import Blueprint, render_template, request, redirect
from src.database.models import (
    get_gwas_trait_name,
    get_locus_info_from_repeat_id,
    check_trait_availability,
    get_available_traits,
    query_manhattan_data,
    get_str_loci_for_trait,
)
from src.plots.manhattan import create_manhattan_plot, create_mini_manhattan_plot
from src.plots.locus import (
    query_allele_data,
    generate_figure_plotly,
    filter_allele_data,
)
from config import Config
import sqlite3
import numpy as np

# Create blueprint
plots_bp = Blueprint("plots", __name__)


@plots_bp.route("/trait_overview/<trait_name>")
def trait_overview(trait_name):
    """Show all Manhattan plots for a trait in a grid layout"""
    # Map trait name for GWAS data lookup
    gwas_trait_name = get_gwas_trait_name(trait_name)

    print(f"Trait overview: {trait_name} -> GWAS: {gwas_trait_name}")

    # Get all STR loci for this trait (use original trait name for locus lookup)
    str_loci = get_str_loci_for_trait(trait_name)

    if not str_loci:
        return (
            render_template(
                "error.html", error=f"No STR loci found for trait '{trait_name}'"
            ),
            404,
        )

    # Check which loci have Manhattan data available (use GWAS trait name)
    available_loci = []
    if os.path.exists(Config.MANHATTAN_DB_PATH):
        try:
            conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)

            for locus in str_loci:
                # Check if this region has data (use GWAS trait name)
                cursor = conn.execute(
                    """
                    SELECT COUNT(*) FROM gwas_variants 
                    WHERE trait_name = ? AND chrom = ? 
                    AND ? BETWEEN (pos - 250000) AND (pos + 250000)
                """,
                    (gwas_trait_name, locus["chrom"], locus["pos"]),
                )

                variant_count = cursor.fetchone()[0]
                if variant_count > 0:
                    locus["variant_count"] = variant_count
                    locus["has_data"] = True
                    available_loci.append(locus)
                else:
                    locus["has_data"] = False

            conn.close()

        except Exception as e:
            print(f"Error checking Manhattan data: {e}")
            # If error, assume no data available
            for locus in str_loci:
                locus["has_data"] = False

    # Generate mini Manhattan plots for each available locus
    plot_data = []
    for locus in available_loci:
        try:
            # Get data for this locus (use GWAS trait name)
            start_pos = max(0, locus["pos"] - 500000)
            end_pos = locus["pos"] + 500000

            data = query_manhattan_data(
                gwas_trait_name, locus["chrom"], start_pos, end_pos
            )

            if not data.empty:
                # Create mini plot
                mini_fig = create_mini_manhattan_plot(
                    data,
                    gwas_trait_name,
                    locus["chrom"],
                    locus["pos"],
                    locus["repeat_id"],
                )

                plot_data.append(
                    {
                        "locus": locus,
                        "plot_json": mini_fig.to_json(),
                        "max_significance": data["neg_log_p"].max()
                        if len(data) > 0
                        else 0,
                        "variant_count": len(data),
                    }
                )

        except Exception as e:
            print(f"Error creating plot for {locus['repeat_id']}: {e}")
            continue

    # Sort by significance (most significant first)
    plot_data.sort(key=lambda x: x["max_significance"], reverse=True)

    return render_template(
        "trait_overview.html",
        trait_name=trait_name,
        plot_data=plot_data,
        total_loci=len(str_loci),
        available_loci=len(available_loci),
    )


@plots_bp.route("/manhattan_plot")
def manhattan_plot_route():
    """Manhattan plot route with improved error handling"""

    # Get parameters
    trait_name = request.args.get("trait")
    repeat_id = request.args.get("repeat_id")
    window_size = request.args.get("window", default=500000, type=int)

    # Map trait name for GWAS data lookup
    gwas_trait_name = get_gwas_trait_name(trait_name)

    print(f"Manhattan plot: {trait_name} -> GWAS: {gwas_trait_name}")

    # Validation
    if not trait_name or not repeat_id:
        return (
            render_template(
                "error.html", error="Missing required parameters: trait and repeat_id"
            ),
            400,
        )

    # Check if databases exist
    if not os.path.exists(Config.LOCUS_DB_PATH):
        return render_template("error.html", error="Locus database not found"), 404

    if not os.path.exists(Config.MANHATTAN_DB_PATH):
        return (
            render_template(
                "error.html",
                error="Manhattan database not found. Please run setup first.",
            ),
            404,
        )

    # Get locus information
    target_chrom, target_pos = get_locus_info_from_repeat_id(repeat_id)
    if target_chrom is None or target_pos is None:
        return (
            render_template(
                "error.html",
                error=f"Could not find locus information for repeat_id {repeat_id}",
            ),
            404,
        )

    # Check trait availability (use GWAS trait name)
    trait_available, trait_status = check_trait_availability(gwas_trait_name)
    if not trait_available:
        available_traits = get_available_traits()
        return (
            render_template(
                "error.html",
                error=f"Trait '{gwas_trait_name}' not available. Status: {trait_status}",
                available_traits=available_traits,
            ),
            404,
        )

    # Query data (use GWAS trait name)
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size

    print(
        f"Querying {gwas_trait_name} data for Chr{target_chrom}:{start_pos}-{end_pos}"
    )
    data = query_manhattan_data(gwas_trait_name, target_chrom, start_pos, end_pos)

    if data.empty:
        return (
            render_template(
                "error.html",
                error=f"No data found for {gwas_trait_name} in region Chr{target_chrom}:{start_pos}-{end_pos}",
            ),
            404,
        )

    # Create plot (use GWAS trait name)
    fig = create_manhattan_plot(
        data, gwas_trait_name, target_chrom, target_pos, window_size
    )

    # Get available traits for dropdown
    available_traits = get_available_traits()

    breadcrumbs = [
        {"name": "Trait Overview", "url": f"/trait_overview/{trait_name}"},
        {"name": f"Chr{target_chrom}:{target_pos:,}", "url": None},
    ]

    # Prepare template data
    template_data = {
        "manhattan_plot_json": fig.to_json(),
        "trait_name": trait_name,
        "repeat_id": repeat_id,
        "locus_chrom": target_chrom,
        "locus_pos": target_pos,
        "window_size": window_size,
        "available_traits": available_traits,
        "data_summary": {
            "total_variants": len(data),
            "min_p": data["p_value"].min() if len(data) > 0 else None,
            "max_p": data["p_value"].max() if len(data) > 0 else None,
            "significant_variants": len(data[data["p_value"] < 5e-8])
            if len(data) > 0
            else 0,
        },
        "breadcrumbs": breadcrumbs,
        "show_trait_overview_link": True,
    }

    return render_template("manhattan_plot.html", **template_data)


# Add this route to your src/routes/plots.py file


@plots_bp.route("/locus_trait_overview/<trait_name>")
def locus_trait_overview(trait_name):
    """Show all locus plots for a trait in a grid layout (similar to trait_overview but for locus plots)"""
    print(f"Locus trait overview: {trait_name}")

    # Get all STR loci for this trait
    str_loci = get_str_loci_for_trait(trait_name)

    if not str_loci:
        return (
            render_template(
                "error.html", error=f"No STR loci found for trait '{trait_name}'"
            ),
            404,
        )

    # Generate mini locus plots for each locus
    plot_data = []
    for locus in str_loci:
        try:
            # Query allele data for this repeat_id
            from src.plots.locus import (
                query_allele_data,
                filter_allele_data,
                create_mini_locus_plot,
            )

            (
                dosage_dict,
                mean_dict,
                ci_dict,
                phenotype,
                db_trait_name,
            ) = query_allele_data(Config.LOCUS_DB_PATH, locus["repeat_id"])

            if dosage_dict and mean_dict:
                # Apply basic filtering
                filtered_dosage, filtered_mean, filtered_ci = filter_allele_data(
                    dosage_dict, mean_dict, ci_dict, count_threshold=50
                )

                if filtered_dosage:
                    # Create mini locus plot
                    mini_fig = create_mini_locus_plot(
                        filtered_dosage,
                        filtered_mean,
                        filtered_ci,
                        trait_name,
                        locus["repeat_id"],
                    )

                    if mini_fig:
                        # Calculate some summary stats
                        allele_count = len(filtered_dosage)
                        total_samples = sum(float(v) for v in filtered_dosage.values())
                        mean_effect = np.mean(
                            [float(v) for v in filtered_mean.values()]
                        )

                        plot_data.append(
                            {
                                "locus": locus,
                                "plot_json": mini_fig.to_json(),
                                "allele_count": allele_count,
                                "total_samples": int(total_samples),
                                "mean_effect": abs(
                                    mean_effect
                                ),  # Use absolute value for sorting
                                "effect_direction": "+" if mean_effect > 0 else "-",
                            }
                        )

        except Exception as e:
            print(f"Error creating locus plot for {locus['repeat_id']}: {e}")
            continue

    # Sort by effect size (largest absolute effect first)
    plot_data.sort(key=lambda x: x["mean_effect"], reverse=True)

    return render_template(
        "locus_trait_overview.html",
        trait_name=trait_name,
        plot_data=plot_data,
        total_loci=len(str_loci),
        available_loci=len(plot_data),
    )


@plots_bp.route("/locus_plot")
def locus_plot_route():
    """Locus plot route for STR allele-phenotype associations with trait switching"""

    # Get parameters
    repeat_id = request.args.get("repeat_id")
    trait_name = request.args.get("trait")
    ci_style = request.args.get("ci_style", default="error_bars")
    color_by_samples = (
        request.args.get("color_by_samples", default="false").lower() == "true"
    )
    count_threshold = request.args.get("count_threshold", default=100, type=int)

    # Validation
    if not repeat_id:
        return (
            render_template(
                "error.html", error="Missing required parameter: repeat_id"
            ),
            400,
        )

    # Check if locus database exists
    if not os.path.exists(Config.LOCUS_DB_PATH):
        return render_template("error.html", error="Locus database not found"), 404

    try:
        # Import the functions from your original locus plotting code
        from src.plots.locus import (
            query_allele_data,
            generate_figure_plotly,
            filter_allele_data,
        )

        # Query allele data
        dosage_dict, mean_dict, ci_dict, phenotype, db_trait_name = query_allele_data(
            Config.LOCUS_DB_PATH, repeat_id
        )

        if dosage_dict is None:
            return (
                render_template(
                    "error.html", error=f"No data found for repeat_id: {repeat_id}"
                ),
                404,
            )

        # Use trait from database if not provided in URL, or use provided trait
        if not trait_name:
            trait_name = db_trait_name or phenotype or "Unknown Trait"

        # Apply filtering
        filtered_dosage, filtered_mean, filtered_ci = filter_allele_data(
            dosage_dict, mean_dict, ci_dict, count_threshold=count_threshold
        )

        if not filtered_dosage:
            return (
                render_template(
                    "error.html",
                    error=f"No data remains after filtering for repeat_id: {repeat_id}",
                ),
                404,
            )

        # Generate the plot
        fig = generate_figure_plotly(
            filtered_dosage,
            filtered_mean,
            filtered_ci,
            trait_name,
            ci_style=ci_style,
            color_by_samples=color_by_samples,
        )

        if fig is None:
            return (
                render_template(
                    "error.html",
                    error=f"Could not generate plot for repeat_id: {repeat_id}",
                ),
                500,
            )

        # Get locus information for display
        target_chrom, target_pos = get_locus_info_from_repeat_id(repeat_id)

        # FIXED: Get ALL available traits instead of just traits for this repeat_id
        # This makes the dropdown work like the Manhattan plot dropdown
        available_traits = get_available_traits()

        # If no traits found, fall back to traits from locus database
        if not available_traits:
            try:
                conn = sqlite3.connect(Config.LOCUS_DB_PATH)
                cursor = conn.execute(
                    """
                    SELECT DISTINCT COALESCE(trait_name, phenotype) as trait
                    FROM locus_data 
                    WHERE (trait_name IS NOT NULL OR phenotype IS NOT NULL)
                    ORDER BY trait
                """
                )
                available_traits = [row[0] for row in cursor.fetchall()]
                conn.close()
            except Exception as e:
                print(f"Error getting available traits: {e}")
                available_traits = [trait_name]  # At least include current trait

        # Map trait name for Manhattan plot link
        gwas_trait_name = get_gwas_trait_name(trait_name)

        # Create breadcrumbs
        breadcrumbs = [
            {"name": "Browse Loci", "url": "/browse_loci"},
            {"name": f"{repeat_id}", "url": None},
        ]

        # Check if Manhattan data is available for this trait
        manhattan_available = False
        if os.path.exists(Config.MANHATTAN_DB_PATH):
            try:
                conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)
                cursor = conn.cursor()
                cursor.execute(
                    "SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?",
                    (gwas_trait_name,),
                )
                variant_count = cursor.fetchone()[0]
                manhattan_available = variant_count > 0
                conn.close()
            except Exception as e:
                print(f"Error checking Manhattan data: {e}")

        # Check if there are other loci for this trait (for trait overview link)
        other_loci_available = False
        try:
            conn = sqlite3.connect(Config.LOCUS_DB_PATH)
            cursor = conn.execute(
                """
                SELECT COUNT(*) FROM locus_data 
                WHERE (trait_name = ? OR phenotype = ?)
                AND repeat_id IS NOT NULL
                AND repeat_id != ?
            """,
                (trait_name, trait_name, repeat_id),
            )
            other_count = cursor.fetchone()[0]
            other_loci_available = other_count > 0
            conn.close()
        except Exception as e:
            print(f"Error checking other loci: {e}")

        # Prepare template data
        template_data = {
            "locus_plot_json": fig.to_json(),
            "repeat_id": repeat_id,
            "trait_name": trait_name,
            "phenotype": phenotype,
            "locus_chrom": target_chrom,
            "locus_pos": target_pos,
            "available_traits": available_traits,  # Now contains ALL available traits
            "ci_style": ci_style,
            "color_by_samples": color_by_samples,
            "count_threshold": count_threshold,
            "breadcrumbs": breadcrumbs,
            "manhattan_available": manhattan_available,
            "gwas_trait_name": gwas_trait_name,
            "other_loci_available": other_loci_available,
            "data_summary": {
                "total_alleles": len(filtered_dosage),
                "original_alleles": len(dosage_dict),
                "filtered_out": len(dosage_dict) - len(filtered_dosage),
            },
        }

        return render_template("locus_plot.html", **template_data)

    except Exception as e:
        return (
            render_template(
                "error.html", error=f"Error generating locus plot: {str(e)}"
            ),
            500,
        )


@plots_bp.route("/test_locus")
def test_locus():
    """Redirect to locus_plot for backwards compatibility"""
    repeat_id = request.args.get("repeat_id")
    trait = request.args.get("trait")

    if repeat_id:
        url_params = f"repeat_id={repeat_id}"
        if trait:
            url_params += f"&trait={trait}"
        return redirect(f"/locus_plot?{url_params}")
    else:
        return redirect("/browse_loci")


@plots_bp.route("/download_manhattan_data")
def download_manhattan_data():
    """Download Manhattan plot data as CSV"""
    from flask import make_response
    import io
    import csv

    # Get parameters
    trait_name = request.args.get("trait")
    repeat_id = request.args.get("repeat_id")
    window_size = request.args.get("window", default=500000, type=int)

    # Map trait name for GWAS data lookup
    gwas_trait_name = get_gwas_trait_name(trait_name)

    # Validation
    if not trait_name or not repeat_id:
        return "Missing required parameters: trait and repeat_id", 400

    # Get locus information
    target_chrom, target_pos = get_locus_info_from_repeat_id(repeat_id)
    if target_chrom is None or target_pos is None:
        return f"Could not find locus information for repeat_id {repeat_id}", 404

    # Query data
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size

    try:
        data = query_manhattan_data(gwas_trait_name, target_chrom, start_pos, end_pos)

        if data.empty:
            return (
                f"No data found for {gwas_trait_name} in region Chr{target_chrom}:{start_pos}-{end_pos}",
                404,
            )

        # Create CSV content
        output = io.StringIO()
        writer = csv.writer(output)

        # Write header
        writer.writerow(
            [
                "chromosome",
                "position",
                "variant_id",
                "p_value",
                "neg_log10_p",
                "beta",
                "se",
            ]
        )

        # Write data
        for _, row in data.iterrows():
            writer.writerow(
                [
                    row["chrom"],
                    row["pos"],
                    row["variant_id"],
                    row["p_value"],
                    row["neg_log_p"],
                    row.get("beta", ""),
                    row.get("se", ""),
                ]
            )

        # Create response
        csv_content = output.getvalue()
        output.close()

        response = make_response(csv_content)
        response.headers["Content-Type"] = "text/csv"
        response.headers[
            "Content-Disposition"
        ] = f"attachment; filename=manhattan_data_{gwas_trait_name}_{repeat_id}_chr{target_chrom}_{target_pos}.csv"

        return response

    except Exception as e:
        return f"Error generating download: {str(e)}", 500
