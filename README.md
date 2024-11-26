## Install

```bash
pip install -e .
```

## Example

```bash
geofetch --mesh "fasting"  --organism human -ds 2023/01/01 --mesh "Diet, Food, and Nutrition" --entry gse --sample rna
```

### Cli

```bash

Usage: geofetch [OPTIONS]

  Fetch GEO data based on user input.

Options:
  -t, --title TEXT                Title(s) of the study or dataset.
  -d, --description TEXT          Description(s) of the study or dataset.
  -o, --organism TEXT             Organism name (default: Homo sapiens).
  -m, --mesh TEXT                 Medical Subject Headings (MeSH) terms.
  -ds, --date-start TEXT          Start date for the data collection (format:
                                  YYYY/MM/DD).
  -de, --date-end TEXT            End date for the data collection (format:
                                  YYYY/MM/DD).
  -e, --entry [gds|gpl|gse|gsm]   Type of entry (gds, gpl, gse, gsm).
  -s, --sample [rna|protein|genomic]
                                  Type of sample (rna, protein, genomic).
  --help                          Show this message and exit.
```
