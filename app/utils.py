from pathlib import Path
import json
import requests
import csv

URL_BASE= "https://pubchem.ncbi.nlm.nih.gov"

def create_folder(folder):
    """Creates a folder."""
    path = Path(folder)
    path.mkdir(parents=True, exist_ok=True)
    return path

def safe_filename(filename):
    """Creates a safe filename"""
    filename = "" if filename is None else str(filename)
    for ch in '<>:/\\|?*':
        filename = filename.replace(ch, "_")
    return filename.strip()


def get_json(url, params=None):
    """GET JSON with 3 retires + timeout. Returns None if failure or timeout"""
    for attempt in range(3):
        try:
            response = requests.get(url, params=params, timeout=30)
            if response.status_code == 200:
                return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Request error for {url}: {e}. Attempt {attempt + 1}/3")   
    return None


def save_json(data, filename, folder):
    """Saves a JSON file in the specified folder."""
    folder_path = create_folder(folder)
    path = folder_path / filename
    path.write_text(json.dumps(data, indent=4), encoding="utf-8") 
    return str(path)

def load_json(filename, folder):
    """Loads a JSON file from the specified folder."""
    folder_path = create_folder(folder)
    path = folder_path / filename

    if not path.exists():
        return None
    
    return json.loads(path.read_text(encoding="utf-8"))


def save_rows_json (rows, filename):
    """Saves a list of rows (dictionaries) as a JSON file."""
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True) 

    with open(path, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=4)

    return str(path)
    
def save_rows_csv (rows, filename):
    """Saves a list of rows (dictionaries) as a CSV file."""
    path = Path(filename) # turns the filename into a Path object
    path.parent.mkdir(parents=True, exist_ok=True) # creates the "parent" folder with the first word in the url
    
    if not rows:
        with open(path, "w", newline="", encoding="utf-8-sig") as f:
            f.write("No data available\n")
            return str(path)
        
    fieldnames = sorted({k for row in rows if isinstance(row,dict) for k in row.keys()})

    with open(path, "w", newline = "", encoding = "utf-8-sig") as f:
        writer = csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            if isinstance(row,dict):
                writer.writerow(row)

    return str(path)

def create_text_file(filename, folder):
    """Creates an empty text file in the specified folder and returns its path."""
    folder_path = create_folder(folder)
    path = folder_path / f"{filename}.txt"

    if not path.exists():
            path.touch() #empty path

    return str(path)

def write_text_file(filename, folder, lines):
    """Opens an existing file and writes on it"""
    folder_path = create_folder(folder)
    path = folder_path / f"{filename}.txt"
    
    text = "\n".join("" if x is None else str(x) for x in lines) + "\n" 
    with open (path, "w") as f: # append to add information to the file without overwritting
        f.write(text)

def read_text_file(filename, folder):
    """Reads a text file and returns its content as a string."""
    folder_path = create_folder(folder)
    path = folder_path / f"{filename}.txt"
    return path.read_text(encoding="utf-8")

def read_file_lines(filename, folder):
    """Reads a text file and returns its content as a list of lines."""
    folder_path = create_folder(folder)
    path = folder_path / f"{filename}.txt"
    return path.read_text(encoding="utf-8").splitlines()