import numpy as np
from Bio import Entrez
import click
import random
import string

def generate_random_email():
    """Generate a random email address."""
    domains = ["example.com", "test.com", "demo.com"]
    username = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
    domain = random.choice(domains)
    return f"{username}@{domain}"


@click.command()
@click.option('-t', '--title', multiple=True, help='Title(s) of the study or dataset.')
@click.option('-d', '--description', multiple=True, help='Description(s) of the study or dataset.')
@click.option('-o', '--organism', default='Homo sapiens', help='Organism name (default: Homo sapiens).')
@click.option('-m', '--mesh', multiple=True, default=['Diabetes Mellitus, Type 2'], help='Medical Subject Headings (MeSH) terms.')
@click.option('-ds', '--date-start', default='2000/01/01', help='Start date for the data collection (format: YYYY/MM/DD).')
@click.option('-de', '--date-end', default='3000', help='End date for the data collection (format: YYYY/MM/DD).')
@click.option('-e', '--entry', default='gds', type=click.Choice(['gds', 'gpl', 'gse', 'gsm']), help='Type of entry (gds, gpl, gse, gsm).')
@click.option('-s', '--sample', default='rna', type=click.Choice(['rna', 'protein', 'genomic']), help='Type of sample (rna, protein, genomic).')

def cli(title, description, organism, mesh, date_start, date_end, sample, entry):
    """Fetch GEO data based on user input."""

    search_term = []


    if sample:
        search_term.append(f'{sample}[Sample Type]')

    if entry:
        search_term.append(f'{entry}[Entry Type]')

    if organism:
        search_term.append(f'{organism}[Organism]')

    mesh_query= []
    if mesh:
        for m in mesh:
            mesh_query.append(f'{m}[MeSH Terms]')
        mesh_query = ' OR '.join(mesh_query)
        search_term.append(f'({mesh_query})')
        
    if date_start and date_end:
        search_term.append(f'({date_start}[Publication Date] : {date_end}[Publication Date])')
        
    description_query = []
    if description:
        for d in description:
            description_query.append(f'{d}[Description]')
        description_query = ' OR '.join(description_query)
        search_term.append(f'({description_query})')
        
    title_query = []
    if title:
        for t in title:
            title_query.append(f'{t}[Title]')
        title_query = ' OR '.join(title_query)
        title_query = f'({title_query})'
        
    search_term = ' AND '.join(search_term)

    Entrez.email = generate_random_email()

    print(f"Using email: {Entrez.email}")
    print(f"Query: {search_term}")
    
    handle = Entrez.esearch(db="gds", term=search_term)
    record = Entrez.read(handle)
    handle.close()
#    print(record)
    print(f"Valid query: {record['QueryTranslation']}")
    print(f"Found {record['Count']} {entry} for above query.")

if __name__ == "__main__":
    cli()
