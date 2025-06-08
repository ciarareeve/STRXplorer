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

# if __name__ == "__main__":
#     app.run(debug=True)

if __name__ == '__main__':
    import os
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)