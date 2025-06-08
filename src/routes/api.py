"""
API route handlers for STRXplorer
"""
import json
from flask import Blueprint, jsonify, current_app
from src.database.models import (
    get_database_stats,
    get_traits_with_loci_data,
    nan_to_null,
)

# Create blueprint
api_bp = Blueprint("api", __name__)


@api_bp.route("/database_status_json")
def database_status_json():
    """Check database status and available traits - RETURNS JSON"""

    status = get_database_stats()

    # Use current_app instead of api_bp for response_class
    return current_app.response_class(
        response=json.dumps(status, default=nan_to_null),
        status=200,
        mimetype="application/json",
    )


@api_bp.route("/api/trait_list")
def api_trait_list():
    """API to get list of traits with available data"""
    try:
        trait_info = get_traits_with_loci_data()

        # Ensure proper data types for JSON serialization
        for trait in trait_info:
            trait["loci_count"] = (
                int(trait["loci_count"]) if trait["loci_count"] is not None else 0
            )
            trait["has_data"] = bool(trait["has_data"])

        return jsonify(trait_info)

    except Exception as e:
        print(f"Error in api_trait_list: {e}")
        return jsonify({"error": str(e)}), 500
