import ollama
import subprocess
import sys
import re

# --- Configuration ---
MODEL_NAME = 'llama4-scout:latest'
WORKER_SCRIPT_NAME = 'generated_fasta_parser.py'
MAX_RETRIES = 5
EXECUTION_TIMEOUT = 300 

BASE_PROMPT = (
    "Write a complete, executable Python script that generates a CRISPRi library.\n\n"
    "Strict Biological & Coding Rules:\n"
    "1. Input Handling: Load 'SpCas9_Indel_2020.csv' containing columns 'sequence' and 'indel'.\n"
    "2. PAM Filtering: You MUST use exactly this syntax: `df_filtered = df[(df['sequence'].str.len() >= 27) & (df['sequence'].str.slice(25, 27).str.upper() == 'GG')].copy()`\n"
    "3. Sequence Counting: Calculate total valid sequences. Print: 'Number of sequences present after NGG filtering: X'. If < 150, exit.\n"
    "4. Value-Based Sampling & Top-Off: You MUST copy and use exactly this code block for sampling:\n"
    "    df_filtered['bin'] = pd.cut(df_filtered['indel'], bins=4, include_lowest=True)\n"
    "    targets = [45, 22, 38, 45] # Bottom, Mid-bottom, Mid-high, Top\n"
    "    sampled_df = pd.DataFrame()\n"
    "    for target, cat in zip(targets, df_filtered['bin'].cat.categories):\n"
    "        bin_df = df_filtered[df_filtered['bin'] == cat]\n"
    "        if len(bin_df) == 0: continue\n"
    "        try:\n"
    "            sampled_df = pd.concat([sampled_df, bin_df.sample(n=target, replace=False)])\n"
    "        except ValueError:\n"
    "            sampled_df = pd.concat([sampled_df, bin_df.sample(n=target, replace=True)])\n"
    "    if len(sampled_df) < 150:\n"
    "        deficit = 150 - len(sampled_df)\n"
    "        sampled_df = pd.concat([sampled_df, df_filtered.sample(n=deficit, replace=False)])\n"
    "5. Sequence Processing: Create a 'guide' column using: `sampled_df['guide'] = sampled_df['sequence'].str.slice(4, 24)`. Transcribe 'T' to 'U' using: `sampled_df['sequence'] = sampled_df['sequence'].str.replace('T', 'U', case=False)`. Do the same for the guide column.\n"
    "6. Output Formatting: Save to 'SpCas9_2020_150_sequence.csv' with columns: sequence, indel, guide.\n\n"
    "Output Constraint: Return ONLY the raw Python code enclosed in a  block. Do not include any markdown explanations."
)

def extract_clean_code(raw_output):
    """Extracts clean Python code from a raw string output by a language model.

    This function attempts to find code enclosed in markdown code blocks (e.g.,
     ...  or  ... ). If no such blocks are found, it
    attempts to clean up the string by removing common markdown prefixes/suffixes.

    Args:
        raw_output: The raw string output from the language model.

    Returns:
        A string containing the extracted and cleaned Python code.
    """
    match = re.search(r'\s*(.*?)\s*', raw_output, re.DOTALL | re.IGNORECASE)
    if match: return match.group(1).strip()
    match2 = re.search(r'\s*(.*?)\s*', raw_output, re.DOTALL)
    if match2: return match2.group(1).strip()
    
    clean = raw_output.strip()
    if clean.lower().startswith(""): clean = clean[9:]
    elif clean.startswith(""): clean = clean[3:]
    if clean.endswith(""): clean = clean[:-3]
    return clean.strip()

def generate_code(prompt):
    """Generates Python code using the Ollama language model.

    Connects to the specified Ollama model, sends the given prompt, and
    processes the model's response to extract clean Python code. Includes
    diagnostic logging for the raw model output.

    Args:
        prompt: The prompt string to send to the language model.

    Returns:
        A string containing the extracted Python code.

    Raises:
        SystemExit: If a fatal error occurs while connecting to Ollama.
    """
    try:
        response = ollama.generate(model=MODEL_NAME, prompt=prompt)
        raw_text = response['response']
        
        # --- CRITICAL DIAGNOSTIC LOGGING ---
        print("\n--- RAW MODEL OUTPUT ---")
        print(raw_text)
        print("------------------------\n")
        
        return extract_clean_code(raw_text)
    except Exception as e:
        print(f"FATAL ERROR connecting to Ollama: {e}")
        sys.exit(1)

def run_worker_script():
    """Executes the worker script and captures its output.

    The script specified by `WORKER_SCRIPT_NAME` is run as a subprocess.
    Its standard output, standard error, and return code are captured.
    A timeout is enforced for script execution.

    Returns:
        A tuple containing:
        - returncode (int): The exit status of the script.
        - stdout (str): The standard output of the script.
        - stderr (str): The standard error of the script.
    """
    try:
        result = subprocess.run([sys.executable, WORKER_SCRIPT_NAME], capture_output=True, text=True, timeout=EXECUTION_TIMEOUT)
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return 1, "", f"TIMEOUT ERROR: Script exceeded {EXECUTION_TIMEOUT} seconds."

def main():
    """Main function for the AI Controller Agent.

    This function orchestrates the code generation and execution process.
    It repeatedly prompts the language model, writes the generated code
    to a file, and executes it. If the script fails, it provides feedback
    to the model and retries, up to `MAX_RETRIES` times.
    """
    print(f"--- Starting AI Controller Agent on {MODEL_NAME} ---")
    current_prompt = BASE_PROMPT
    
    for attempt in range(1, MAX_RETRIES + 1):
        print(f"\n[Attempt {attempt}/{MAX_RETRIES}] Generating code from model...")
        code = generate_code(current_prompt)
        
        # --- STRICT LENGTH VALIDATION ---
        if not code or len(code) < 10:
            print("\n[!] The extracted code is empty or too short. Forcing a retry...")
            stderr = "You failed to output a valid Python script. You must output raw code."
            returncode = 1
            stdout = ""
        else:
            with open(WORKER_SCRIPT_NAME, "w") as f:
                f.write(code)
            print(f"Code saved to {WORKER_SCRIPT_NAME}. Executing...")
            returncode, stdout, stderr = run_worker_script()
        
        if stdout:
            print(f"\n--- Worker Standard Output ---\n{stdout.strip()}")
            
        if stderr or returncode != 0:
            print(f"\n--- Worker Encountered an Issue ---\n{stderr.strip()}")
            if attempt == MAX_RETRIES:
                print("\n[!] Max retries reached.")
                break
            
            # Keep the biological rules in the context during retries
            current_prompt = (
                "ORIGINAL RULES:\n" + BASE_PROMPT + "\n\n"
                "YOUR PREVIOUS OUTPUT FAILED WITH ERROR:\n" + stderr + "\n\n"
                "Rewrite the script to fix this. Return ONLY the code."
            )
        else:
            print(f"\n[SUCCESS] Script executed cleanly on attempt {attempt}.")
            break

if __name__ == "__main__":
    main()
