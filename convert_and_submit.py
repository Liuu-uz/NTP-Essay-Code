import os
import time
import requests

# Define directories
base_dir = os.path.abspath("results")
fasta_files = [os.path.join(base_dir, f) for f in os.listdir(base_dir) if f.endswith(".fasta")]

# URL for form submission (action URL from the form)
url = "http://biomine.cs.vcu.edu/servers/biomine.php?name=NsitePred"

# Define your email address
your_email = "liuzhijing09@gmail.com"  # Replace this with your actual email address

# Function to submit the form
def submit_fasta(file_path, email):
    with open(file_path, 'r') as f:
        fasta_content = f.read()

    # Form data to simulate form submission
    form_data = {
        "seq": fasta_content,  # Input field for sequence
        "email1": email,       # Input field for email
    }

    # Send the POST request
    response = requests.post(url, data=form_data, allow_redirects=True)

    # Check the response
    if response.status_code == 200:
        print(f"Submission successful for {file_path}")
        # if "Your submission is now being processed" in response.text:
        print(f"Submission successful for {file_path}. Submission is being processed.")
        # else:
        #     print(f"Submission for {file_path} succeeded, but processing confirmation not detected.")
    elif response.status_code in [301, 302]:
        redirect_url = response.headers.get("Location")
        print(f"Submission successful for {file_path}. Redirected to: {redirect_url}")
    else:
        print(f"Failed to process {file_path}. Status code: {response.status_code}")
        print(f"Response headers: {response.headers}")
        print(f"Response content: {response.text[:500]}")  # Print the first 500 characters for debugging

# Process each FASTA file
for fasta_file in fasta_files:
    print(f"Processing {fasta_file}...")
    submit_fasta(fasta_file, your_email)
    time.sleep(7)  # Wait for 7 seconds between each submission