from pathlib import Path
import json
import requests
import csv

URL_base = "https://pubchem.ncbi.nlm.nih.gov"

# Check problems and timeout
def get_json(url):
    for attempt in range(3): # Try 3 times to retrieve the data
        print(f"\nAttempting to retrieve data from URL: {url} (Attempt {attempt + 1}/3)")
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                print(f"Successfully retrieved data from URL in attempt {attempt + 1}")
                print("Printing JSON...")
                return response.json()
            else:
                print(f"Request failed with status code {response.status_code}. Attempt {attempt + 1}/3")
        # Catch any request exceptions
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}. Attempt {attempt + 1}/3")   
    print("Failed to retrieve data after 3 attempts.")
    return None


# Save JSON to file GENERAL USE
def save_json(json_file, filename, folder):
    Path(folder).mkdir(parents=True, exist_ok=True)
    path = Path(folder) / filename
    path.write_text(json.dumps(json_file, indent=4)) # writes the JSON data to the file with indentation for readability
    return str(path)

def make_subfolder(parent_folder, child_folder):
    path = Path(parent_folder) / child_folder
    path.mkdir(parents=True, exist_ok=True)
    return str(path)

# Find if the file exists and return it
def load_json(filename, folder):
    Path(folder).mkdir(parents=True, exist_ok=True) # creates the folder if it doesn't exist
    path = Path(folder) / filename
    if path.exists():
        print(f"File {path} already exists. Loading from file.")
        print("\n")
        return json.loads(path.read_text())
    else:
        # Give the option to the user to click download files?? THINK ABOUT IT
        # TODO
        print(f"File {filename} does not exist. Please run the code to create it.")
        return None
    


def save_rows_json (rows, compound, filename):
    path = Path(filename) # turns the filename into a Path object
    path.parent.mkdir(parents=True, exist_ok=True) # creates the "parent" folder with the first word in the url
    with open(path, "w") as f:
        json.dump(rows, f, indent=4)
        return str(path)
    

def save_rows_csv (rows, compound, filename):
    path = Path(filename) # turns the filename into a Path object
    path.parent.mkdir(parents=True, exist_ok=True) # creates the "parent" folder with the first word in the url

    file_path = path

    if not rows:
        print("No rows to save in CSV format for file: "+filename)
        with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
            f.write("No data available\n")
            return str(file_path)
        
    fieldnames = sorted({k for row in rows for k in row.keys()})
    with open(file_path, "w", newline = "", encoding = "utf-8-sig") as f:
        writer = csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return str(file_path)

# Save JSON to file GENERAL USE
def save_json(json_file, filename, folder):
    Path(folder).mkdir(parents=True, exist_ok=True)
    path = Path(folder) / filename
    path.write_text(json.dumps(json_file, indent=4)) # writes the JSON data to the file with indentation for readability
    return str(path)

# Create a file
def create_file(filename, folder):
    Path(folder).mkdir(parents=True, exist_ok=True)
    path = Path(folder) / f"{filename}.txt"

    if path.exists():
            print("File already exists, loading from cache...")
            return str(path)

    path.touch() # Creates an empty path
    return str(path)

# Open an existing file and write on it
def write_file(filename, folder, lines):
    path = Path(folder) / f"{filename}.txt"
    text = "\n".join("" if x is None else str(x) for x in lines) + "\n" # TODO cambiar esto??
    with open (path, "w") as f: # append to add information to the file without overwritting
        f.write(text)

# Reads a txt file and returns the line read
def read_file(filename, folder):
    path = Path(folder) / f"{filename}.txt"
    with open (path, "r") as f:
        line = f.readline()
    return line

