import os, json, urllib.request, time, glob

api_key = os.environ.get("GOOGLE_API_KEY")
if not api_key:
    print("Error: GOOGLE_API_KEY environment variable not set.")
    exit(1)

url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-3.1-pro-preview:generateContent?key={api_key}"

for filepath in glob.glob("*.py"):
    if filepath == "api_docstring.py" or filepath.endswith("_standard.py"):
        continue
        
    print(f"Processing {filepath}...")
    with open(filepath, 'r') as f:
        code = f.read()

    prompt = "Add Google-style docstrings. Do NOT change logic. Output ONLY raw python code, no markdown.\n\n" + code
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    req = urllib.request.Request(url, data=json.dumps(payload).encode('utf-8'), headers={'Content-Type': 'application/json'})

    try:
        with urllib.request.urlopen(req) as res:
            data = json.loads(res.read().decode())
            text = data['candidates'][0]['content']['parts'][0]['text']
            
            # Clean up markdown safely
            text = text.replace("```python", "").replace("```", "").strip()
            
            with open(filepath.replace(".py", "_standard.py"), "w") as f:
                f.write(text + "\n")
    except Exception as e:
        print(f"Error on {filepath}: {e}")

    time.sleep(10)
