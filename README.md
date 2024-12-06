## Install

```bash
pip install -e .
```

## Example

```bash
geofetch --entry gse --sample rna --sample mpss  --sample sage  --organism human  --mesh 'diet' --mesh 'body weight' --mesh-operator AND
```

### Cli

```bash
Usage: geofetch [OPTIONS]

  Query GEO database for series, samples, and datasets. version 4.3.0

Options:
  --count / --print-records       print counts or print records.  [default:
                                  print-records]
  -o, --organism TEXT             Organism name (default: Homo sapiens).
                                  [default: Homo sapiens]
  -ds, --date-start TEXT          Start date for the data collection (format:
                                  YYYY/MM/DD).  [default: 2000/01/01]
  -de, --date-end TEXT            End date for the data collection (format:
                                  YYYY/MM/DD).  [default: 3000]
  -e, --entry [gse|gsm]           Entry type to search. GPL and GDS support
                                  might be added if/when needed.  [default:
                                  gsm]
  -m, --mesh TEXT                 Medical Subject Headings (MeSH) terms.
                                  [default: Diabetes Mellitus, Type 2]
  -mo, --mesh-operator [OR|AND]   Operator for Medical Subject Headings (MeSH)
                                  terms.  [default: OR]
  -s, --sample [rna|mpss|sage|protein|genomic|any]
                                  Type of sample.  [default: rna]
  -t, --title TEXT                Title(s) of the study or dataset.
  -d, --description TEXT          Description(s) of the study or dataset.
  --log-level [WARNING|INFO|DEBUG|ERROR]
                                  Logging level in terms of urgency  [default:
                                  WARNING]
  -fw, --file-write               flag to enable to write to file
  -fn, --file-name TEXT           Output file name.  [default:
                                  geofetch_20241206]
  -ft, --file-type [csv|excel|parquet]
                                  Output file type.  [default: csv]
  --help                          Show this message and exit.
