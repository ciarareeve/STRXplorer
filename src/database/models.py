"""
Database models and query functions for STRXplorer
"""
import sqlite3
import pandas as pd
import os
import math
from typing import Optional, Tuple, List
from config import Config


def get_gwas_trait_name(trait_name):
    """
    Map trait names from locus database to Manhattan database format
    Handles the naming differences between the two databases
    """
    # Mapping from locus DB names (with 'mean_') to Manhattan DB names (without 'mean_')
    trait_mapping = {
        'mean_corpuscular_haemoglobin': 'corpuscular_haemoglobin',
        'mean_corpuscular_volume': 'corpuscular_volume', 
        'mean_platelet_volume': 'platelet_volume',
        'mean_sphered_cell_volume': 'sphered_cell_volume',
        # Add any other mappings as needed
        # 'locus_db_name': 'manhattan_db_name'
    }
    
    # Return mapped name if it exists, otherwise return original name
    return trait_mapping.get(trait_name, trait_name)


def get_locus_info_from_repeat_id(
    repeat_id: str,
) -> Tuple[Optional[str], Optional[int]]:
    """Get chromosome and position for a repeat_id"""
    try:
        conn = sqlite3.connect(Config.LOCUS_DB_PATH)
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT chrom, pos
            FROM locus_data
            WHERE repeat_id = ?
        """,
            (repeat_id,),
        )

        result = cursor.fetchone()
        conn.close()

        if result:
            chrom, pos = result
            # Ensure chromosome format is consistent (remove 'chr' prefix if present)
            chrom_clean = str(chrom).replace("chr", "")
            return chrom_clean, int(pos)
        else:
            return None, None

    except Exception as e:
        print(f"Error getting locus info: {e}")
        return None, None


def query_manhattan_data(
    trait_name: str, chrom: str, start_pos: int, end_pos: int
) -> pd.DataFrame:
    """Query Manhattan plot data from database"""
    if not os.path.exists(Config.MANHATTAN_DB_PATH):
        return pd.DataFrame()

    try:
        conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)

        query = """
            SELECT chrom, pos, variant_id, p_value, neg_log_p, beta, se
            FROM gwas_variants
            WHERE trait_name = ? AND chrom = ? AND pos BETWEEN ? AND ?
            AND p_value IS NOT NULL AND neg_log_p IS NOT NULL
            ORDER BY pos
        """

        df = pd.read_sql_query(
            query, conn, params=(trait_name, chrom, start_pos, end_pos)
        )
        conn.close()

        return df

    except Exception as e:
        print(f"Error querying Manhattan data: {e}")
        return pd.DataFrame()


def check_trait_availability(trait_name: str) -> Tuple[bool, str]:
    """Check if a trait is available in the database"""
    if not os.path.exists(Config.MANHATTAN_DB_PATH):
        return False, "Manhattan database not found"

    try:
        conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT total_variants FROM trait_metadata 
            WHERE trait_name = ?
        """,
            (trait_name,),
        )

        result = cursor.fetchone()
        conn.close()

        if result:
            return True, f"Available ({result[0]:,} variants)"
        else:
            return False, "Not in database"

    except Exception as e:
        return False, f"Error: {e}"


def get_available_traits() -> List[str]:
    """Get list of traits available in the database"""
    if not os.path.exists(Config.MANHATTAN_DB_PATH):
        return []

    try:
        conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)
        cursor = conn.cursor()

        cursor.execute("SELECT trait_name FROM trait_metadata ORDER BY trait_name")
        traits = [row[0] for row in cursor.fetchall()]

        conn.close()
        return traits

    except Exception as e:
        print(f"Error getting available traits: {e}")
        return []


def get_traits_with_loci_data() -> List[dict]:
    """Get traits with data availability from locus database"""
    try:
        conn = sqlite3.connect(Config.LOCUS_DB_PATH)
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

            if os.path.exists(Config.MANHATTAN_DB_PATH):
                manhattan_conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)
                manhattan_cursor = manhattan_conn.execute(
                    "SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?",
                    (gwas_trait_name,),
                )
                variant_count = manhattan_cursor.fetchone()[0]
                has_manhattan_data = variant_count > 0
                manhattan_conn.close()

            trait_info.append(
                {
                    "trait": trait,
                    "gwas_trait": gwas_trait_name,
                    "loci_count": loci_count,
                    "has_data": has_manhattan_data,
                    "display_name": trait.replace("_", " ").title(),
                    "status": "Available" if has_manhattan_data else "No GWAS Data",
                }
            )

        conn.close()
        return trait_info

    except Exception as e:
        print(f"Error loading traits: {e}")
        return []


def get_all_str_loci() -> List[dict]:
    """Get all STR loci with their trait associations"""
    try:
        conn = sqlite3.connect(Config.LOCUS_DB_PATH)
        cursor = conn.execute(
            """
            SELECT repeat_id, chrom, pos, motif, ref_len,
                   COALESCE(trait_name, phenotype) as trait
            FROM locus_data 
            WHERE repeat_id IS NOT NULL
            AND (trait_name IS NOT NULL OR phenotype IS NOT NULL)
            ORDER BY chrom, pos
        """
        )

        loci_info = []
        for row in cursor.fetchall():
            repeat_id, chrom, pos, motif, ref_len, trait = row

            # Handle None values properly
            loci_info.append(
                {
                    "repeat_id": repeat_id,
                    "chrom": str(chrom).replace("chr", "") if chrom else "Unknown",
                    "pos": int(pos) if pos is not None else 0,
                    "motif": motif if motif else "Unknown",
                    "ref_len": float(ref_len) if ref_len is not None else 0.0,
                    "trait": trait if trait else "Unknown",
                    "display_name": trait.replace("_", " ").title()
                    if trait
                    else "Unknown",
                    "location": f"Chr{str(chrom).replace('chr', '')}:{pos:,}"
                    if chrom and pos
                    else "Unknown",
                }
            )

        conn.close()
        return loci_info

    except Exception as e:
        print(f"Error loading STR loci: {e}")
        return []


def get_str_loci_for_trait(trait_name: str) -> List[dict]:
    """Get all STR loci for a specific trait"""
    try:
        conn = sqlite3.connect(Config.LOCUS_DB_PATH)
        cursor = conn.execute(
            """
            SELECT repeat_id, chrom, pos, motif, ref_len 
            FROM locus_data 
            WHERE (trait_name = ? OR phenotype = ?)
            AND repeat_id IS NOT NULL
            ORDER BY chrom, pos
        """,
            (trait_name, trait_name),
        )

        str_loci = []
        for row in cursor.fetchall():
            repeat_id, chrom, pos, motif, ref_len = row
            str_loci.append(
                {
                    "repeat_id": repeat_id,
                    "chrom": str(chrom).replace("chr", "") if chrom else "Unknown",
                    "pos": int(pos) if pos is not None else 0,
                    "motif": motif if motif else "Unknown",
                    "ref_len": float(ref_len) if ref_len is not None else 0.0,
                    "location": f"chr{str(chrom).replace('chr', '')}:{pos:,}"
                    if chrom and pos
                    else "Unknown",
                }
            )
        conn.close()
        return str_loci

    except Exception as e:
        print(f"Error accessing locus database: {e}")
        return []


def get_database_stats() -> dict:
    """Get database statistics"""
    stats = {
        "locus_db_exists": os.path.exists(Config.LOCUS_DB_PATH),
        "manhattan_db_exists": os.path.exists(Config.MANHATTAN_DB_PATH),
        "available_traits": get_available_traits()
        if os.path.exists(Config.MANHATTAN_DB_PATH)
        else [],
        "manhattan_db_path": Config.MANHATTAN_DB_PATH,
        "locus_db_path": Config.LOCUS_DB_PATH,
        "total_variants": 0,
    }

    # Get database stats if available
    if stats["manhattan_db_exists"]:
        try:
            conn = sqlite3.connect(Config.MANHATTAN_DB_PATH)
            cursor = conn.cursor()

            # Total variants
            cursor.execute("SELECT COUNT(*) FROM gwas_variants")
            stats["total_variants"] = cursor.fetchone()[0]

            # Trait metadata
            cursor.execute(
                """
                SELECT trait_name, total_variants, min_p_value, max_p_value
                FROM trait_metadata
                ORDER BY trait_name
            """
            )

            trait_stats = []
            for row in cursor.fetchall():
                trait_name, total_variants, min_p_value, max_p_value = row

                # Convert NaN values to None (which becomes null in JSON)
                trait_stats.append(
                    {
                        "trait_name": trait_name,
                        "total_variants": total_variants
                        if total_variants is not None
                        else 0,
                        "min_p_value": min_p_value
                        if min_p_value is not None and not math.isnan(min_p_value)
                        else None,
                        "max_p_value": max_p_value
                        if max_p_value is not None and not math.isnan(max_p_value)
                        else None,
                    }
                )

            stats["trait_stats"] = trait_stats
            conn.close()

        except Exception as e:
            stats["error"] = str(e)

    return stats


def nan_to_null(obj):
    """Convert NaN values to None for JSON serialization"""
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return obj
