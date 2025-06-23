"""
STRXplorer - Main Flask Application
Modular architecture with blueprints for better organization
"""
from flask import Flask
from config import config
from src.routes.main import main_bp
from src.routes.plots import plots_bp
from src.routes.api import api_bp


def create_app(config_name="default"):
    """
    Application factory pattern

    Args:
        config_name: Configuration to use ('development', 'production', or 'default')

    Returns:
        Flask application instance
    """
    app = Flask(__name__)

    # Load configuration
    app.config.from_object(config[config_name])

    # Register blueprints
    app.register_blueprint(main_bp)
    app.register_blueprint(plots_bp)
    app.register_blueprint(api_bp)

    return app


# Create the app instance
app = create_app()

# Add this right after creating the app
import os
print(f"DEBUG: Current working directory: {os.getcwd()}")
print(f"DEBUG: Files in current directory: {os.listdir('.')}")
print(f"DEBUG: Files in data directory: {os.listdir('data') if os.path.exists('data') else 'data directory not found'}")
print(f"DEBUG: locus_data.db exists: {os.path.exists('data/locus_data.db')}")
print(f"DEBUG: manhattan_data.db exists: {os.path.exists('data/manhattan_data.db')}")

# if __name__ == "__main__":
#     app.run(debug=True)

if __name__ == '__main__':
    import os
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)