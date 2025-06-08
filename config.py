"""
Configuration settings for STRXplorer application
"""
import os


class Config:
    """Base configuration class"""

    # Database paths (updated to use data/ directory)
    LOCUS_DB_PATH = "data/locus_data.db"
    MANHATTAN_DB_PATH = "data/manhattan_data.db"

    # Flask settings
    SECRET_KEY = os.environ.get("SECRET_KEY") or "str-xplorer-secret-key"

    # Trait name mapping for GWAS data
    TRAIT_NAME_MAPPING = {
        "corpuscular_haemoglobin": "mean_corpuscular_haemoglobin",
        "corpuscular_volume": "mean_corpuscular_volume",
        "platelet_volume": "mean_platelet_volume",
        "sphered_cell_volume": "mean_sphered_cell_volume",
    }


class DevelopmentConfig(Config):
    """Development configuration"""

    DEBUG = True


class ProductionConfig(Config):
    """Production configuration"""

    DEBUG = False


# Default configuration
config = {
    "development": DevelopmentConfig,
    "production": ProductionConfig,
    "default": DevelopmentConfig,
}
