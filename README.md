# Snakemake workflow: Kingfisher_workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/trickovicmatija/kingfisher_workflow/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for downloading whole projects together with their metadata using [kingfisher-download](https://github.com/wwood/kingfisher-download). It extends its functionality by making possible to download samples from the same BioProject in parallel. Additionally, it merges runs coming from one sample together, to make further analysis easier.


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

you need to install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).


Run as
```

snakemake -d output_dir -j1 separate_samples --config BioProject=PRJNA363003 --use-conda --resources ncbi_connection=4 

```

Run as
```

snakemake -d output_dir -j10 --config BioProject=PRJNA363003 --use-conda --resources ncbi_connections=4

```


# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.
