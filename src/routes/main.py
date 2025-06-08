"""
Main route handlers for STRXplorer
"""
from flask import Blueprint, render_template
from src.database.models import (
    get_available_traits,
    get_database_stats,
    get_traits_with_loci_data,
    get_all_str_loci,
    get_gwas_trait_name,
)
import os
import sqlite3


# Create blueprint
main_bp = Blueprint("main", __name__)


@main_bp.route("/")
def home():
    """Main landing page"""
    # Get available traits (same as Manhattan page)
    available_traits = get_available_traits()

    # Get database stats
    db_stats = get_database_stats()
    stats = {
        "available_traits": available_traits,
        "total_variants": db_stats.get("total_variants", 0),
        "total_traits": len(available_traits),
        "loaded_traits": len(available_traits),
    }

    return render_template("home.html", available_traits=available_traits, stats=stats)


# Update your browse_traits route in src/routes/main.py


@main_bp.route("/browse_traits")
def browse_traits():
    """Browse all available traits with enhanced plot type checking"""
    try:
        # Get traits with STR loci
        conn = sqlite3.connect("data/locus_data.db")
        cursor = conn.execute(
            """
            SELECT DISTINCT 
                COALESCE(trait_name, phenotype) as trait,
                COUNT(repeat_id) as loci_count
            FROM locus_data 
            WHERE repeat_id IS NOT NULL
            AND (trait_name IS NOT NULL OR phenotype IS NOT NULL)
            GROUP BY trait
            HAVING loci_count > 0
            ORDER BY loci_count DESC, trait
        """
        )

        trait_info = []
        for row in cursor.fetchall():
            trait, loci_count = row

            # Map trait name and check if trait has Manhattan data
            gwas_trait_name = get_gwas_trait_name(trait)
            has_manhattan_data = False

            if os.path.exists("data/manhattan_data.db"):
                manhattan_conn = sqlite3.connect("data/manhattan_data.db")
                manhattan_cursor = manhattan_conn.execute(
                    "SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?",
                    (gwas_trait_name,),
                )
                variant_count = manhattan_cursor.fetchone()[0]
                has_manhattan_data = variant_count > 0
                manhattan_conn.close()

            # Check if trait has locus plot data (JSON format with allele data)
            has_locus_data = False
            try:
                locus_cursor = conn.execute(
                    """
                    SELECT data_json FROM locus_data 
                    WHERE (trait_name = ? OR phenotype = ?)
                    AND repeat_id IS NOT NULL 
                    AND data_json IS NOT NULL
                    LIMIT 1
                """,
                    (trait, trait),
                )

                sample_row = locus_cursor.fetchone()
                if sample_row:
                    import json

                    data = json.loads(sample_row[0])
                    # Check if it has the expected locus data structure
                    has_locus_data = "sample_count_per_summed_length" in data or any(
                        key.startswith("mean_") for key in data.keys()
                    )
            except Exception as e:
                print(f"Error checking locus data for {trait}: {e}")

            trait_info.append(
                {
                    "trait": trait,
                    "gwas_trait": gwas_trait_name,
                    "loci_count": loci_count,
                    "has_data": has_manhattan_data,  # For backward compatibility
                    "has_manhattan_data": has_manhattan_data,
                    "has_locus_data": has_locus_data,
                    "display_name": trait.replace("_", " ").title(),
                    "status": "Available" if has_manhattan_data else "No GWAS Data",
                }
            )

        conn.close()
        return render_template("browse_traits.html", traits=trait_info)

    except Exception as e:
        return render_template("error.html", error=f"Error loading traits: {e}"), 500


@main_bp.route("/browse_loci")
def browse_loci():
    """Browse STR loci with search/filter capabilities"""
    try:
        loci_info = get_all_str_loci()

        # Group by chromosome for better organization
        loci_by_chrom = {}
        for locus in loci_info:
            chrom = locus["chrom"]
            if chrom not in loci_by_chrom:
                loci_by_chrom[chrom] = []
            loci_by_chrom[chrom].append(locus)

        # Get unique traits for the filter dropdown
        unique_traits = list(
            set([locus["trait"] for locus in loci_info if locus["trait"] != "Unknown"])
        )

        # ADD THESE MISSING VARIABLES:

        # 1. Create trait_mapping (maps internal trait names to GWAS trait names)
        trait_mapping = {}
        for trait in unique_traits:
            trait_mapping[trait] = get_gwas_trait_name(trait)

        # 2. Get available Manhattan traits from the database
        available_manhattan_traits = set()
        if os.path.exists("data/manhattan_data.db"):
            try:
                manhattan_conn = sqlite3.connect("data/manhattan_data.db")
                cursor = manhattan_conn.execute(
                    "SELECT DISTINCT trait_name FROM gwas_variants"
                )
                available_manhattan_traits = set(row[0] for row in cursor.fetchall())
                manhattan_conn.close()
            except Exception as e:
                print(f"Error loading Manhattan traits: {e}")
                available_manhattan_traits = set()

        return render_template(
            "browse_loci.html",
            loci_by_chrom=loci_by_chrom,
            total_loci=len(loci_info),
            unique_traits=sorted(unique_traits),
            trait_mapping=trait_mapping,  # ADD THIS
            available_manhattan_traits=available_manhattan_traits,
        )  # ADD THIS

    except Exception as e:
        print(f"Error loading STR loci: {e}")  # Add this for debugging
        return render_template("error.html", error=f"Error loading STR loci: {e}"), 500


@main_bp.route("/about")
def about():
    """About page describing the project and team"""
    return render_template("about.html")


@main_bp.route("/database_status")
def database_status():
    """Database status with GUI - RETURNS HTML"""
    return render_template("database_status.html")


@main_bp.route("/setup_instructions")
def setup_instructions():
    """Show setup instructions"""
    return render_template("setup_instructions.html")
