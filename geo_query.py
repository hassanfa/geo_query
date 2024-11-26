import numpy as np
from Bio import Entrez
import click
import random
import string

def fetch_geo_data(tissue, organism, disease, email, start_date=None, end_date=None):
    # Set the email for Entrez
    Entrez.email = email
    
    # Construct the search term
    search_term = f"{tissue}[Tissue] AND {organism}[Organism] AND {disease}[Disease]"
    
    # Add date range if provided
    if start_date and end_date:
        search_term += f" AND ({start_date}[Publication Date] : {end_date}[Publication Date])"
    
    # Search GEO
    handle = Entrez.esearch(db="gds", term=search_term, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    
    # Fetch GEO data using the IDs obtained from the search
    geo_ids = record["IdList"]
    list_of_GEO = []
    
    for geo_id in geo_ids:
        handle = Entrez.efetch(db="gds", id=geo_id, retmode="xml")
        geo_data = Entrez.read(handle)
        handle.close()
        list_of_GEO.append(geo_data)
    
    return list_of_GEO

def generate_random_email():
    """Generate a random email address."""
    domains = ["example.com", "test.com", "demo.com"]
    username = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
    domain = random.choice(domains)
    return f"{username}@{domain}"

@click.command()
@click.option('--tissue', prompt='Tissue type', help='Type of tissue to search for.')
@click.option('--organism', prompt='Organism', help='Organism to search for.')
@click.option('--disease', prompt='Disease', help='Disease to search for.')
@click.option('--email', default=generate_random_email(), help='Your email for NCBI Entrez.')
@click.option('--start-date', default=None, help='Start date (YYYY/MM/DD) for publication date filter.')
@click.option('--end-date', default=None, help='End date (YYYY/MM/DD) for publication date filter.')
def main(tissue, organism, disease, email, start_date, end_date):
    """Fetch GEO data based on user input."""
    print(f"Using email: {email}")  # Display the email being used
    data = fetch_geo_data(tissue, organism, disease, email, start_date, end_date)
    
    # Print the results
    for entry in data:
        print(entry)

if __name__ == "__main__":
    main()
