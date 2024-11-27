import numpy as np
from Bio import Entrez
import click
import random
import string


def generate_random_email():
    """Generate a random email address."""
    domains = ["example.com", "test.com", "demo.com"]
    username = ''.join(
        random.choices(string.ascii_lowercase + string.digits, k=8))
    domain = random.choice(domains)
    return f"{username}@{domain}"


def build_query(terms, query_type, operator='OR'):
    query = []
    operator = f' {operator} '
    if terms:
        for term in terms:
            query.append(f'{term}[{query_type}]')
        return f'({operator.join(query)})'
    return ''


@click.command()
@click.option('-o',
              '--organism',
              default='Homo sapiens',
              help='Organism name (default: Homo sapiens).')
@click.option('-ds',
              '--date-start',
              default='2000/01/01',
              help='Start date for the data collection (format: YYYY/MM/DD).')
@click.option('-de',
              '--date-end',
              default='3000',
              help='End date for the data collection (format: YYYY/MM/DD).')
@click.option('-e',
              '--entry',
              default='gse',
              type=click.Choice(['gds', 'gpl', 'gse', 'gsm']),
              help='Type of entry (gds, gpl, gse, gsm).')
@click.option('-m',
              '--mesh',
              multiple=True,
              default=['Diabetes Mellitus, Type 2'],
              help='Medical Subject Headings (MeSH) terms.')
@click.option('-s',
              '--sample',
              multiple=True,
              default=['rna'],
              type=click.Choice(['rna', 'mpss', 'sage', 'protein', 'genomic']),
              help='Type of sample.')
@click.option('-t',
              '--title',
              multiple=True,
              help='Title(s) of the study or dataset.')
@click.option('-d',
              '--description',
              multiple=True,
              help='Description(s) of the study or dataset.')
def cli(title, description, organism, mesh, date_start, date_end, sample,
        entry):
    """Fetch GEO data based on user input."""

    search_term = []

    search_term.append(f'{entry}[Entry Type]')
    search_term.append(f'{organism}[Organism]')
    search_term.append(
        f'({date_start}[Publication Date] : {date_end}[Publication Date])')
    search_term.append(
        build_query(terms=sample, query_type='Sample Type', operator='OR'))
    search_term.append(
        build_query(terms=mesh, query_type='MeSH Terms', operator='OR'))
    if description:
        search_term.append(
            build_query(terms=description,
                        query_type='Description',
                        operator='OR'))
    if title:
        search_term.append(
            build_query(terms=title, query_type='Title', operator='OR'))

    search_term = ' AND '.join(search_term)

    Entrez.email = generate_random_email()

    print(f"Using email: {Entrez.email}")
    print(f"Query: {search_term}")


#    handle = Entrez.esearch(db="gds", term=search_term)
#    record = Entrez.read(handle)
#    handle.close()
#    print(record)
#    print(f"Valid query: {record['QueryTranslation']}")
#    print(f"Found {record['Count']} {entry} for above query.")

if __name__ == "__main__":
    cli()
