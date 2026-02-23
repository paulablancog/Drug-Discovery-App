from pathlib import Path
import json
import requests

URL_base = "https://pubchem.ncbi.nlm.nih.gov"

# Check problems and timeout
def get_json(url):
    for attempt in range(3): # Try 3 times to retrieve the data
        print(f"Attempting to retrieve data from URL: {url} (Attempt {attempt + 1}/3)")
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                print(f"Successfully retrieved data from URL in attempt {attempt + 1}")
                print("Printing index JSON...")
                print("\n")
                print("Top keys: ", list(response.json().keys()))
                print("Record keys: ", list(response.json()["Record"].keys()))
                print("Number of top sections: ", len(response.json()["Record"].get("Section", [])))
                return response.json()
            else:
                print(f"Request failed with status code {response.status_code}. Attempt {attempt + 1}/3")
        # Catch any request exceptions
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}. Attempt {attempt + 1}/3")   
    print("Failed to retrieve data after 3 attempts.")
    return None


# Save index JSON to file
def save_json(index_json, filename):
    path = Path(filename) # turns the filename into a Path object
    path.parent.mkdir(parents=True, exist_ok=True) # creates the "parent" folder with the first word in the url
    path.write_text(json.dumps(index_json, indent=4)) # writes the JSON data to the file with indentation for readability
    return str(path)


# Find if the file exists and return it
def load_index_json(filename, folder="indexes"):
    Path(folder).mkdir(parents=True, exist_ok=True) # creates the folder if it doesn't exist
    file = Path(folder) / filename
    if file.exists():
        print(f"File {filename} already exists. Loading from file.")
        return json.loads(file.read_text())
    else:
        # Give the option to the user to click download files?? THINK ABOUT IT
        # TODO
        print(f"File {filename} does not exist. Please run the code to create it.")
        return None