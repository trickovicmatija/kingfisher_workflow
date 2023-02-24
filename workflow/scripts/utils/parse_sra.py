# from ..color_logger import logger
import logging

logger = logging.getLogger(__file__)
import pandas as pd


Expected_library_values = {
    "library_selection": "RANDOM",
    "library_strategy": "WGS",
    "library_source": "METAGENOMIC",
}


def load_and_validate_runinfo_table(path):

    RunTable = pd.read_csv(path, sep="\t", index_col=0)

    # validate sra table
    format_error = False

    # check if all headers are present
    Expected_headers = [
        "read1_length_average",
        "read2_length_average",
        "library_source",
        "library_selection",
        "library_strategy",
        "sample_accession",
    ]
    for header in Expected_headers:
        if not header in RunTable.columns:

            logger.error(f"Didn't found expected header {header}")
            format_error = True

    if not all(RunTable.index.str[1:3] == "RR"):
        logger.error("Expect runs as index, e.g. [E,S,D]RR000")
        format_error = True

    if not (RunTable.sample_accession.str[1:3]=="RS").all():
        logger.error("sample_accession should be something like [E,S]RS' ")
        format_error = True


    if format_error:
        logger.error("RunTable {} is not valid. Abort.".format(path))
        exit(1)

    return RunTable


def filter_runinfo(RunTable, ignore_paired=False):

    logger.info(
        f"Start with {RunTable.shape[0]} runs from {RunTable.sample_accession.unique().shape[0]} samples"
    )

    # Filter out reads that are not metagenomics

    for key in ["library_source"]:

        Nruns_before = RunTable.shape[0]
        All_values = RunTable[key].unique()
        RunTable = RunTable.loc[RunTable[key].str.lower().str.contains(Expected_library_values[key].lower())]

        Difference = Nruns_before - RunTable.shape[0]

        if Difference > 0:

            logger.info(
                f"Runs have the folowing values for {key}: {', '.join(All_values)}\n"
                f"Select only runs that contain {Expected_library_values[key]} in {key}, "
                f"Filtered out {Difference} runs"
            )

    for key in ["library_selection", "library_strategy"]:

        Nruns_before = RunTable.shape[0]
        All_values = RunTable[key].unique()
        if any(RunTable[key] != Expected_library_values[key]):

            logger.warning(
                f"Runs have the folowing values for {key}: {', '.join(All_values)}\n"
                f"Usually I expect {key} == {Expected_library_values[key]} "
            )

   

    # Final
    if RunTable.shape[0] > 0:
        logger.info(
            f"Selected {RunTable.shape[0]} runs from {RunTable.sample_accession.unique().shape[0]} samples"
        )

    else:
        logger.critical("No runs left after filtering. Abort.")
        exit(1)

    return RunTable


def validate_merging_runinfo(path):

    RunTable = load_and_validate_runinfo_table(path)

    # If each run is from a different sample, merging is not necessary
    if RunTable.shape[0] == RunTable.sample_accession.unique().shape[0]:
        return RunTable

    # Cannot merge if different platforms
    problematic_samples = []
    for sample, df in RunTable.groupby("sample_accession"):
        if not all(df.Platform == df.Platform.iloc[0]):
            problematic_samples.append(sample)

    if len(problematic_samples) > 0:
        logger.error(
            f"You attemt to merge runs from the same sample. "
            f"But for {len(problematic_samples)} samples the runs are sequenced with different platforms and should't be merged.\n"
            f"Please resolve the the abiguity in the table {path} and rerun the command.\n"
        )

        exit(1)

    # Warn if samples are not identical for the follwing columns
    Expected_same_values = ["experiment_accession", "model", "library_name"]
    for key in Expected_same_values:

        problematic_samples = []
        for sample, df in RunTable.groupby("sample_accession"):
            if not all(df[key] == df[key].iloc[0]):
                problematic_samples.append(sample)

        if len(problematic_samples) > 0:
            if len(problematic_samples) > 5:
                problematic_samples_list = " ".join(problematic_samples[:3] + ["..."])
            else:
                problematic_samples_list = " ".join(problematic_samples)

                logger.warn(
                    "You attemt to merge runs from the same sample. "
                    f"But for {len(problematic_samples)} samples the runs have different {key}: {problematic_samples_list}\n"
                    f"You can modify the table {path} and rerun the command.\n"
                )

    logger.info("I will automatically merge runs from the same sample.")

    return RunTable
