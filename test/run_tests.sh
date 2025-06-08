#!/bin/bash
# Automated testing script for STRXplorer

echo "🚀 STRXplorer Automated Testing"
echo "================================"

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 is not installed or not in PATH"
    exit 1
fi

# Check if required packages are installed
echo "📦 Checking dependencies..."
python3 -c "import requests, flask, sqlite3, pandas, plotly" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "❌ Missing required Python packages"
    echo "Please install requirements:"
    echo "   pip install -r requirements.txt"
    exit 1
fi

# Set Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Check if databases exist
if [ ! -f "data/locus_data.db" ] || [ ! -f "data/manhattan_data.db" ]; then
    echo "⚠️  Warning: Database files not found in data/ directory"
    echo "   Expected: data/locus_data.db and data/manhattan_data.db"
    echo "   Tests may fail without proper data"
fi

# Start Flask server in background
echo "🌐 Starting Flask server..."
python3 app.py &
SERVER_PID=$!

# Wait for server to start
echo "⏳ Waiting for server to start..."
sleep 3

# Check if server started successfully
if ! kill -0 $SERVER_PID 2>/dev/null; then
    echo "❌ Failed to start Flask server"
    exit 1
fi

# Run tests
echo "Running endpoint tests..."
python3 test_endpoints.py

# Capture test exit code
TEST_EXIT_CODE=$?

# Stop the server
echo "Stopping Flask server..."
kill $SERVER_PID
wait $SERVER_PID 2>/dev/null

# Report final status
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "✅ All tests completed successfully!"
else
    echo "❌ Some tests failed (exit code: $TEST_EXIT_CODE)"
fi

exit $TEST_EXIT_CODE