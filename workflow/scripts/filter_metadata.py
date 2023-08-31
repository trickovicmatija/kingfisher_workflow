import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

### SCRIPT START ###


import pandas as pd

from utils.parse_sra import load_and_validate_runinfo_table, filter_runinfo

# Parse runtable
RunTable = load_and_validate_runinfo_table(snakemake.input[0])

# Filter runtable
RunTable_filtered = filter_runinfo(RunTable, ignore_paired=False)


RunTable_filtered.to_csv(snakemake.output[0], sep="\t")