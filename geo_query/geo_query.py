"""
Main function to query GEO using a click interface
"""

import os
import random
import string
import logging
from datetime import datetime

import click
import polars as pl
from Bio import Entrez
from xlsxwriter import Workbook

from geo_query.version import __version__ as version


def write_df_to_file(df, filename, filetype, logger, search_term, worksheet="geofetch"):
    """
    Write a Polars DataFrame to a file.
    """
    comments = [
        f"# geofetch version: {version}",
        f"# date: {datetime.now().strftime('%Y%m%d')}",
        f"# query: {search_term}",
    ]
    try:
        if filetype.lower() == "csv":
            with open(filename, mode="a") as f:
                for line in comments:
                    f.write(f"{line}\n")
                df.write_csv(f)
        else:
            with Workbook(filename) as wb:
                # comments start at the start of the position
                sheet_position = [0, 0]
                df_comments = pl.DataFrame(comments)
                df_comments.write_excel(
                    workbook=wb,
                    worksheet=worksheet,
                    position=tuple(sheet_position),
                    include_header=False,
                )

                # update sheet position to rows after comments
                sheet_position[0] += df_comments.height
                df.write_excel(
                    workbook=wb,
                    worksheet=worksheet,
                    position=tuple(sheet_position),
                    include_header=False,
                )

        logger.info(f"File successfully written to {filename}")
    except Exception as e:
        logger.error(f"Error writing file: {str(e)}")
        raise e


def default_filename():
    """
    creates a default filename with a today date suffix
    """
    return f"geofetch_{datetime.now().strftime('%Y%m%d')}"


def validate_filename(ctx, param, value):
    """
    validates filename and handles the errors e.g., file exists
    """
    if value is None:
        value = default_filename()

    if os.path.exists(value):
        # pass error to context to handle in the main
        ctx.filename_error = f"Error: The file '{value}' already exists."
        return value

    return value


def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange.
    """

    def document(func):
        func.__doc__ = docstring
        return func

    return document


class EntrezGDS:
    """
    A class to interface Entrez to query GEO and process the results
    """

    def __init__(self, email, logger):
        Entrez.email = email
        self.logger = logger

    def esearch(self, search_term, retmax):
        """
        Query GEO database for a search term and retmax numbers
        """
        handle = Entrez.esearch(db="gds", term=search_term, retmax=retmax)
        self.logger.debug(f"{handle}")
        record = Entrez.read(handle)
        self.logger.debug(record)
        self.logger.debug(f"Valid query: {record['QueryTranslation']}")
        handle.close()
        return record

    def esummary(self, search_id):
        """
        Query GEO for the IDs that have been found previously using esearch
        """
        handle = Entrez.esummary(db="gds", id=search_id)
        self.logger.debug(f"{handle}")
        record = Entrez.read(handle)
        handle.close()
        return record

    def process_record(self, search_term, mesh, count):
        """
        Query GEO using the esearch and search_term for count number of items
        and return a polars dataframe
        """
        record = self.esearch(search_term=search_term, retmax=count)
        if int(count) > 9999:
            self.logger.warning(
                f"{count} results requested. Only the first 10 000 will be returned\n"
                f"(GEO query limitations). Limit your search terms"
            )

        if len(record["IdList"]) == 0:
            self.logger.warning("Search result returned zero output")
            return None

        summary = self.esummary(search_id=record["IdList"])
        entry_type = summary[1]["entryType"].lower()

        if entry_type == "gsm":
            return self._process_gsm(summary, mesh)
        return self._process_gse(summary, mesh)

    def _process_gsm(self, summary, mesh):
        df = pl.DataFrame(
            [
                {
                    "GSE": s["GSE"],
                    "GPL": s["GPL"],
                    "GSMtaxon": s["taxon"],
                    "GSM": s["Accession"],
                }
                for s in summary
            ]
        )

        df = self._add_mesh_column(df, mesh)
        df = self._explode_columns(df, ["GSE", "mesh"])
        df = self._prefix_column(df, "GSE")

        return df

    def _process_gse(self, summary, mesh):
        df = pl.DataFrame(
            [
                {
                    "GSE": s["Accession"],
                    "GPL": s["GPL"],
                    "GSEtaxon": s["taxon"],
                    "GSM": [sample["Accession"] for sample in s["Samples"]],
                    "SampleCount": len(s["Samples"]),
                }
                for s in summary
            ]
        )

        df = self._add_mesh_column(df, mesh)
        df = self._explode_columns(df, ["GSE", "mesh"])

        return df

    def _add_mesh_column(self, df, mesh):
        return df.with_columns(pl.lit(";".join(mesh)).alias("mesh"))

    def _explode_columns(self, df, columns):
        for col in columns:
            df = df.with_columns([pl.col(col).str.split(";")]).explode(col)
        return df

    def _prefix_column(self, df, column):
        return df.with_columns((pl.lit(column) + pl.col(column)).alias(column))


def generate_random_email():
    """
    A generator for random email address to use for Entrez
    """
    domains = ["example.com", "test.com", "demo.com"]
    username = "".join(random.choices(string.ascii_lowercase + string.digits, k=8))
    domain = random.choice(domains)
    return f"{username}@{domain}"


def build_query(terms, query_type, operator="OR"):
    """
    A query builder given list of terms, query type, and aggregation operator
    """
    query = []
    operator = f" {operator} "
    if terms:
        for term in terms:
            query.append(f"{term}[{query_type}]")
        return f"({operator.join(query)})"
    return ""


def initialize_logger(log_level):
    """
    Initalizes logger and encapsulates it
    """
    logging.basicConfig(
        level=logging.ERROR,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    return logger


@click.command()
@click.option(
    "--count/--print-records",
    "count",
    default=False,
    is_flag=True,
    show_default=True,
    help="print counts or print records.",
)
@click.option(
    "-o",
    "--organism",
    default=["Homo sapiens"],
    multiple=True,
    show_default=True,
    help="Organism name (default: Homo sapiens).",
)
@click.option(
    "-ds",
    "--date-start",
    default="2000/01/01",
    show_default=True,
    help="Start date for the data collection (format: YYYY/MM/DD).",
)
@click.option(
    "-de",
    "--date-end",
    default="3000",
    show_default=True,
    help="End date for the data collection (format: YYYY/MM/DD).",
)
@click.option(
    "-e",
    "--entry",
    default="gsm",
    show_default=True,
    type=click.Choice(["gse", "gsm"]),
    help="Entry type to search. GPL and GDS support might be added if/when needed.",
)
@click.option(
    "-m",
    "--mesh",
    multiple=True,
    show_default=True,
    default=["Diabetes Mellitus, Type 2"],
    help="Medical Subject Headings (MeSH) terms.",
)
@click.option(
    "-mo",
    "--mesh-operator",
    default="OR",
    show_default=True,
    type=click.Choice(["OR", "AND"]),
    help="Operator for Medical Subject Headings (MeSH) terms.",
)
@click.option(
    "-s",
    "--sample",
    multiple=True,
    default=["rna"],
    show_default=True,
    type=click.Choice(["rna", "mpss", "sage", "protein", "genomic"]),
    help="Type of sample.",
)
@click.option("-t", "--title", multiple=True, help="Title(s) of the study or dataset.")
@click.option(
    "-d", "--description", multiple=True, help="Description(s) of the study or dataset."
)
@click.option(
    "--log-level",
    default="WARNING",
    type=click.Choice(["WARNING", "INFO", "DEBUG", "ERROR"], case_sensitive=False),
    help="Logging level in terms of urgency",
    show_default=True,
)
@click.option(
    "-fw",
    "--file-write",
    is_flag=True,
    default=False,
    help="flag to enable to write to file",
)
@click.option(
    "-fn",
    "--file-name",
    default=default_filename(),
    show_default=True,
    callback=validate_filename,
    help="Output file name.",
)
@click.option(
    "-ft",
    "--file-type",
    default="csv",
    type=click.Choice(["csv", "excel"]),
    show_default=True,
    help="Output file type.",
)
@add_doc(f"Query GEO database for series, samples, and datasets. version {version}")
@click.pass_context
def cli(
    ctx,
    title,
    description,
    organism,
    mesh,
    mesh_operator,
    date_start,
    date_end,
    sample,
    entry,
    log_level,
    count,
    file_write,
    file_name,
    file_type,
):
    """Fetch GEO data based on user input."""
    logger = initialize_logger(getattr(logging, log_level))

    # if output file exists, throw an error and exit
    if hasattr(ctx, "filename_error"):
        logger.error(ctx.filename_error)
        ctx.exit(1)
    else:
        logger.info(f"Output will be written to {file_name}")

    search_term = []

    search_term.append(f"{entry}[Entry Type]")

    search_term.append(
        build_query(terms=organism, query_type="Organism", operator="OR")
    )

    search_term.append(
        f"({date_start}[Publication Date] : {date_end}[Publication Date])"
    )
    search_term.append(
        build_query(terms=sample, query_type="Sample Type", operator="OR")
    )
    search_term.append(
        build_query(terms=mesh, query_type="MeSH Terms", operator=mesh_operator)
    )
    if description:
        search_term.append(
            build_query(terms=description, query_type="Description", operator="OR")
        )
    if title:
        search_term.append(build_query(terms=title, query_type="Title", operator="OR"))

    search_term = " AND ".join(search_term)

    logger.debug(f"Query: {search_term}")

    entrez_gds = EntrezGDS(email=generate_random_email(), logger=logger)

    record = entrez_gds.esearch(search_term=search_term, retmax=0)

    if count:
        click.echo(f"{record['Count']}")
    else:

        df = entrez_gds.process_record(
            search_term=search_term, mesh=mesh, count=record["Count"]
        )

        if file_write:
            write_df_to_file(
                df=df,
                search_term=search_term,
                filename=file_name,
                filetype=file_type,
                logger=logger,
            )

        if not df is None:
            pl.Config.set_tbl_rows(-1)
            click.echo(df.shape)
            click.echo(df.head())
            pl.Config.set_tbl_rows(10)


if __name__ == "__main__":
    cli()
