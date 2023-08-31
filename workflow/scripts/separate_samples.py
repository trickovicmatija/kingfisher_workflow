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

from pathlib import Path
from utils import parse_sra


runinfo_file = snakemake.input[0]
output_folder = Path(snakemake.output[0])

# The script assumes that all are either paired end or single end
# TODO better handle paired /single end samples.
# e.g. drop singletons if other library is paired

RunTable = parse_sra.validate_merging_runinfo(runinfo_file)


# All biosamples
Samples = RunTable.sample_accession.unique()



# check if all are paired end



if 'read2_length_average' not in RunTable.columns:
    paired=False
else:

    samples_with_missing_r2 = RunTable.index[RunTable.read2_length_average.isnull()]
    if len(samples_with_missing_r2) == RunTable.shape[0]:
        paired = False
    elif len(samples_with_missing_r2)==0:
        paired = True
    else: 

        logger.error(
            f"Your library layout is not consistent, please check your runtable {runinfo_file}"
        )
        exit(1)




output_folder.mkdir()


for sample in RunTable["sample_accession"].unique():
    df = RunTable.query("sample_accession == @sample")
    
    
    with open(output_folder / (sample+".txt") ,"w") as outf:
        outf.write("\n".join( df.index ) +'\n')
    



"""


# check if sample.tsv already exist

sample_table_file = working_dir / "samples.tsv"
# delete samples.tsv if it exists and overwrite is set
if sample_table_file.exists():
        #sample_table_file.unlink()

    logger.error(
        f"{sample_table_file} already exists, I dare not to overwrite it. "
    )
    exit(1)

# create sample table

sample_table = pd.DataFrame(index=Samples)

SRA_READ_PATH = SRA_subfolder.relative_to(working_dir) / "Samples"

if not paired:
    sample_table["R1"] = sample_table.index.map(
        lambda s: str(SRA_READ_PATH / f"{s}/{s}.fastq.gz")
    )
else:

    sample_table["R1"] = sample_table.index.map(
        lambda s: str(SRA_READ_PATH / f"{s}/{s}_1.fastq.gz")
    )
    sample_table["R2"] = sample_table.index.map(
        lambda s: str(SRA_READ_PATH / f"{s}/{s}_2.fastq.gz")
    )

prepare_sample_table_for_atlas(
    sample_table, reads_are_QC=skip_qc, outfile=str(sample_table_file)
)

logger.info(f"Prepared sample table with {sample_table.shape[0]} samples")



"""
