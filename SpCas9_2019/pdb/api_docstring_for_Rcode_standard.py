"""Processes R scripts in the current directory to add roxygen2 documentation.

This script searches for all `.R` and `.r` files in the working directory,
excluding those that already have a `_standard` suffix. It sends the content
of each file to the Google Gemini API to generate standard roxygen2 documentation
headers without altering the original R logic. The documented code is then
saved to a new file with the `_standard` suffix. A 10-second delay is added
between requests to respect rate limits.

Environment Variables:
    GOOGLE_API_KEY: The Google Gemini API key used for authentication.

Raises:
    SystemExit: If the GOOGLE_API_KEY environment variable is not set.
"""
import os, json, urllib.request, time, glob

api_key = os.environ.get("GOOGLE_API_KEY")
if not api_key:
    print("Error: GOOGLE_API_KEY environment variable not set.")
    exit(1)

url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={api_key}"

# 2. Target both .R and .r files
for filepath in glob.glob("*.[rR]"):
    if filepath.endswith("_standard.R") or filepath.endswith("_standard.r"):
        continue
        
    print(f"Processing {filepath}...")
    with open(filepath, 'r') as f:
        code = f.read()

    # 3. Analytical prompt specifically requesting roxygen2 standard
    prompt = (
        "Add standard roxygen2 documentation headers (using #') to the following R code. "
        "Do NOT change any logic, variable names, or arguments. "
        "Output ONLY raw R code, no markdown explanations.\n\n"
    ) + code
    
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    req = urllib.request.Request(url, data=json.dumps(payload).encode('utf-8'), headers={'Content-Type': 'application/json'})

    try:
        with urllib.request.urlopen(req) as res:
            data = json.loads(res.read().decode())
            text = data['candidates'][0]['content']['parts'][0]['text']
            
            # 4. Clean up R-specific markdown tags
            text = text.replace("R", "").replace("r", "").replace("", "").strip()
            
            # Save the file preserving the original extension casing
            out_name = filepath[:-2] + "_standard" + filepath[-2:]
            with open(out_name, "w") as f:
                f.write(text + "\n")
                
    except Exception as e:
        print(f"Error on {filepath}: {e}")

    # 10-second delay to respect the free tier rate limit
    time.sleep(10)
