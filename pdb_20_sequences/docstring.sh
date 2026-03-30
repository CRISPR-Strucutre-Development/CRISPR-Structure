# Set your API Key as an environment variable (copy it from your Google AI Studio)
export GOOGLE_API_KEY="AIzaSyCTxVV4XrVsU_gCvuBLY7H6L0Mp3n6GxBc"

for file in *.py; do
    echo "Processing $file..."
    # 1. Use the current Gemini 3 Flash model
    # 2. Use -p for non-interactive automation
    gemini -m gemini-3-flash-preview -p "Add Google-style docstrings. Do NOT change logic." < "$file" > "${file%.py}_standard.py"
    
    # Stay under free-tier limits
    sleep 10
done
