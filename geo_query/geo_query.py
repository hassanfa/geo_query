import click
import random
import string
import logging
import polars as pl

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


class EntrezGDS:

    def __init__(self, email, logger):
        Entrez.email = email
        self.logger = logger

    def esearch(self, search_term, retmax):
        handle = Entrez.esearch(db="gds", term=search_term, retmax=retmax)
        self.logger.debug(f"{handle}")
        record = Entrez.read(handle)
        handle.close()
        return record

    def esummary(self, search_id):
        handle = Entrez.esummary(db="gds", id=search_id)
        self.logger.debug(f"{handle}")
        record = Entrez.read(handle)
        handle.close()
        return record


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
@click.option('--count/--print-records',
              'count',
              default=False,
              is_flag=True,
              show_default=True,
              help='print counts or print records.')
@click.option('-o',
              '--organism',
              default=['Homo sapiens'],
              multiple=True,
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
@click.option(
    '-e',
    '--entry',
    default='gsm',
    show_default=True,
    type=click.Choice(['gse', 'gsm']),
    help=
    'Entry type to search. GPL and GDS support might be added if/when needed.')
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

    search_term.append(
        build_query(terms=organism, query_type='Organism', operator='OR'))

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

    logger.debug(f"Query: {search_term}")

    entrez_gds = EntrezGDS(email=generate_random_email(), logger=logger)

    record = entrez_gds.esearch(search_term=search_term, retmax=0)
    logger.debug(record)
    logger.debug(f"Valid query: {record['QueryTranslation']}")
    if count:
        click.echo(f"{record['Count']}")
    else:
        max_print = 20
        if int(record["Count"]) > max_print:
            logger.info(
                f"There are {record['Count']} items found. Printing only the first {max_print}"
            )
        record = entrez_gds.esearch(search_term=search_term,
                                    retmax=record["Count"][0:max_print - 1])

        click.echo(f"all samples: {record['IdList'][0:49]}")

        summary = entrez_gds.esummary(search_id=record['IdList'])

        # check if summary actually has an output!
        if summary[1]['entryType'].lower() == 'gsm':
            df = pl.DataFrame([{
                "GSE": s["GSE"],
                "GPL": s["GPL"],
                "GSMtaxon": s["taxon"],
                "GSM": s["Accession"]
            } for s in summary])

            df = df.with_columns(pl.lit(";".join(mesh)).alias('mesh'))

            for col in ['GSE', 'mesh']:
                df = df.with_columns([pl.col(col).str.split(';')]).explode(col)

            for col in ['GSE']:
                df = df.with_columns((pl.lit(col) + pl.col(col)).alias(col))
        elif summary[1]['entryType'].lower() == 'gse':
            df = pl.DataFrame([{
                "GSE":
                s["Accession"],
                "GPL":
                s["GPL"],
                "GSEtaxon":
                s["taxon"],
                "GSM": [sample['Accession'] for sample in s['Samples']]
            } for s in summary])

            df = df.with_columns(pl.lit(";".join(mesh)).alias('mesh'))

            for col in ['GSE', 'mesh']:
                df = df.with_columns([pl.col(col).str.split(';')]).explode(col)

        pl.Config.set_tbl_rows(-1)
        click.echo(df)
        pl.Config.set_tbl_rows(10)


if __name__ == "__main__":
    cli()
