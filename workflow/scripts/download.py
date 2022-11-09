import glob,os
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
output_dir = f"{snakemake.config['output_dir']}/reads"

def run_kingfisher(sample_list, output_dir, threads, log):
    """Run kingfisher on a list of samples"""
    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Run kingfisher
    shell(
        f"cd {output_dir} && kingfisher get --run-identifiers-list {sample_list} --download-threads {threads}  -m ena-ftp -f fastq.gz &> {log}"
    )


with open(snakemake.input.sample, 'r') as sample_file:
    runs = [run.rstrip("\n") for run in sample_file.readlines()]

if len(runs) == 1:
    run_kingfisher(snakemake.input.sample, output_dir, snakemake.threads, snakemake.log[0])
    if len(glob.glob(f"{output_dir}/{runs[0]}*.fastq.gz")) == 2:
        shell(f"mv {output_dir}/{runs[0]}_1.fastq.gz {output_dir}/{snakemake.wildcards.sample}_1.fastq.gz")
        shell(f"mv {output_dir}/{runs[0]}_2.fastq.gz {output_dir}/{snakemake.wildcards.sample}_2.fastq.gz")
    elif len(glob.glob(f"{output_dir}/{runs[0]}*.fastq.gz")) == 1:
        shell(f"mv {output_dir}/{runs[0]}.fastq.gz {output_dir}/{snakemake.wildcards.sample}.fastq.gz")

    shell(f"touch {snakemake.output[0]}")
elif len(runs) > 1:
    run_kingfisher(snakemake.input.sample, f"{snakemake.config['tmp_dir']}/{snakemake.wildcards.sample}", snakemake.threads, snakemake.log[0])
    shell(f"cat {snakemake.config['tmp_dir']}/{snakemake.wildcards.sample}/*_1.fastq.gz > {output_dir}/{snakemake.wildcards.sample}_1.fastq.gz")
    shell(f"cat {snakemake.config['tmp_dir']}/{snakemake.wildcards.sample}/*_2.fastq.gz > {output_dir}/{snakemake.wildcards.sample}_2.fastq.gz")
    shell(f"touch {snakemake.output[0]}")
else:
    logging.error("No runs found for this sample")