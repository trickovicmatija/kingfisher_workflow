import glob,os
import pandas as pd
from snakemake.shell import shell
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

# Read in the sample sheet
run_accessions = ("SRR","ERR","DRR")
biosample_accesions = ("SAMEA","SAMN","SAMD","SAMG","SAMR","SAMS","SAMT","SAMX")

os.makedirs(snakemake.output.samples, exist_ok=True)

sample_sheet = pd.read_csv(snakemake.input.metadata, sep="\t")

sample_column = ""
run_column = ""

if snakemake.config['filter_samples']:
    allowed_values = tuple(snakemake.config['allowed_values'])
    removed_runs = sample_sheet[~sample_sheet['library_selection'].str.contains("|".join(allowed_values),case=False)].iloc[:,0].values
    if len(removed_runs) > 0:
        sample_sheet = sample_sheet[sample_sheet['library_selection'].str.contains("|".join(allowed_values),case=False)]
        logging.warning(f"Removed {len(removed_runs)} runs from the sample sheet.")
        logging.warning(f"Removed runs: {removed_runs}")
    else:
        logging.info("No runs removed from the sample sheet.")


for column in sample_sheet.columns:

    potential_biosample = sample_sheet[column].astype('str').str.startswith(biosample_accesions)
    potential_run = sample_sheet[column].astype('str').str.startswith(run_accessions)

    if potential_biosample.sum() == potential_biosample.shape[0]:
        sample_column = column

    if potential_run.sum() == potential_run.shape[0]:
        run_column = column
    
    if sample_column and run_column:
        break

for sample in sample_sheet[sample_column].unique():
    df = sample_sheet[sample_sheet[sample_column] == sample]
    df[[run_column]].to_csv(f"{snakemake.output.samples}/{sample}.txt", sep="\t", index=False,header=False)