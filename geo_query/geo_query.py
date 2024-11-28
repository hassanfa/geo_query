import click
import random
import string
import logging

from Bio import Entrez
from enum import StrEnum
from typing import List


def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange.
    """

    def document(func):
        func.__doc__ = docstring
        return func

    return document


def generate_random_email():
    """
    A generator for random email address to use for Entrez
    """
    domains = ["example.com", "test.com", "demo.com"]
    username = ''.join(
        random.choices(string.ascii_lowercase + string.digits, k=8))
    domain = random.choice(domains)
    return f"{username}@{domain}"


def build_query(terms, query_type, operator='OR'):
    """
    A query builder given list of terms, query type, and aggregation operator
    """
    query = []
    operator = f' {operator} '
    if terms:
        for term in terms:
            query.append(f'{term}[{query_type}]')
        return f'({operator.join(query)})'
    return ''


@click.command()
@click.option('-C',
              '--count',
              default=False,
              is_flag=True,
              show_default=True,
              help='Print only count')
@click.option('-o',
              '--organism',
              default='Homo sapiens',
              show_default=True,
              help='Organism name (default: Homo sapiens).')
@click.option('-ds',
              '--date-start',
              default='2000/01/01',
              show_default=True,
              help='Start date for the data collection (format: YYYY/MM/DD).')
@click.option('-de',
              '--date-end',
              default='3000',
              show_default=True,
              help='End date for the data collection (format: YYYY/MM/DD).')
@click.option('-e',
              '--entry',
              default='gse',
              show_default=True,
              type=click.Choice(['gds', 'gpl', 'gse', 'gsm']),
              help='Type of entry (gds, gpl, gse, gsm).')
@click.option('-m',
              '--mesh',
              multiple=True,
              show_default=True,
              default=['Diabetes Mellitus, Type 2'],
              help='Medical Subject Headings (MeSH) terms.')
@click.option('-mo',
              '--mesh-operator',
              default='OR',
              show_default=True,
              type=click.Choice(['OR', 'AND']),
              help='Operator for Medical Subject Headings (MeSH) terms.')
@click.option('-s',
              '--sample',
              multiple=True,
              default=['rna'],
              show_default=True,
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
@click.option(
    "--log-level",
    default='INFO',
    type=click.Choice(['INFO', 'DEBUG']),
    help="Logging level in terms of urgency",
    show_default=True,
)
@add_doc("Query GEO database for series, samples, and datasets.}")
def cli(title, description, organism, mesh, mesh_operator, date_start,
        date_end, sample, entry, log_level, count):
    """Fetch GEO data based on user input."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    search_term = []

    search_term.append(f'{entry}[Entry Type]')
    search_term.append(f'{organism}[Organism]')
    search_term.append(
        f'({date_start}[Publication Date] : {date_end}[Publication Date])')
    search_term.append(
        build_query(terms=sample, query_type='Sample Type', operator='OR'))
    search_term.append(
        build_query(terms=mesh,
                    query_type='MeSH Terms',
                    operator=mesh_operator))
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

    logger.debug(f"Using email: {Entrez.email}")
    logger.debug(f"Query: {search_term}")

    handle = Entrez.esearch(db="gds", term=search_term)
    logger.debug(f"{handle}")
    record = Entrez.read(handle)
    handle.close()
    logger.debug(record)
    logger.debug(f"Valid query: {record['QueryTranslation']}")
    if count:
        click.echo(f"{record['Count']}")
    else:
        click.echo(f"Found {record['Count']} {entry} for above query.")


if __name__ == "__main__":
    cli()
