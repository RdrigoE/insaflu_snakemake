from typing import Literal, TypedDict

from validation.validation import Validator


class ProjectSettings(TypedDict):
    only_samples: bool
    abricate: bool
    project_name: str
    fasta_reference: str
    gb_reference: str
    illumina_consensus: Literal["snippy"] | Literal["iVar"]
    primers_fasta: str | None
    get_minor_variants: bool


class ProjectPath(TypedDict):
    project_dir: str
    reference_dir: str
    primer_dir: str | None


def validate_project_settings(project_settings: ProjectSettings, paths: ProjectPath):
    errors = {}

    if not isinstance(project_settings["only_samples"], bool):
        errors["only_samples"] = "You must provide true or false to only_samples."
        return errors

    if not isinstance(project_settings["abricate"], bool):
        errors["abricate"] = "You must provide true or false to abricate."

    if project_settings["only_samples"] is True:
        return errors

    if not Validator.string(project_settings["project_name"]):
        errors[
            "project_name"
        ] = "You choose to run in project mode, please provide a name without spaces to the field 'project_name'."

    if not Validator.is_fasta(
        paths["reference_dir"] + project_settings["fasta_reference"]
    ):
        errors["fasta_reference"] = "Please provide a valid FASTA file."

    if not Validator.is_genbank(
        paths["reference_dir"] + project_settings["gb_reference"]
    ):
        errors["gb_reference"] = "Please provide a valid GenBank file."

    if not Validator.same_identifiers(
        paths["reference_dir"] + project_settings["fasta_reference"],
        paths["reference_dir"] + project_settings["gb_reference"],
    ):
        errors[
            "identifiers"
        ] = "Please provide a FASTA and GenBank with the same identifiers."

    if not Validator.includes(
        project_settings["illumina_consensus"], ["snippy", "iVar"]
    ):
        errors[
            "illumina_consensus"
        ] = f"You provided the tool '{project_settings['illumina_consensus']}' but currently the pipeline only supports 'snippy' and 'iVar'.\nChange the illumina_consensus field."

    if project_settings["primers_fasta"] is not None and not Validator.is_fasta(
        paths["primer_dir"] + project_settings["primers_fasta"]
    ):
        errors["primers"] = "Please provide valid primers in the FASTA format."

    if project_settings[
        "primers_fasta"
    ] is not None and not Validator.check_file_in_the_directory(
        paths["primer_dir"], project_settings["primers_fasta"] +
            ".pair_information.tsv"
    ):
        errors[
            "primers_files"
        ] = "Please provide the pair_information primer file with the same name as the primer with the extension '.pair_information.tsv'."

    if not isinstance(project_settings["get_minor_variants"], bool):
        errors[
            "get_minor_variants"
        ] = "You must provide true or false to get_minor_variants."

    return errors

    # REFERENCE_GB = f"../user/references/{user_metadata['gb_reference']}"
    # assert is_genbank(REFERENCE_GB), "Reference genbank file is not in genbank format"
    # REFERENCE_FASTA = f"../user/references/{user_metadata['fasta_reference']}"
    # assert is_fasta(REFERENCE_FASTA), "Reference fasta file is not in fasta format"
    # assert same_identifiers(
    #     REFERENCE_FASTA, REFERENCE_GB
    # ), "Reference fasta and genbank files do not have the same identifiers"
    # REFERENCE_NAME = re.findall("(?<=references/)(.*?)(?=.fasta)", REFERENCE_FASTA)[0]
    #
    #
    #
    #
    # PROJECT_NAME = user_metadata["project_name"]
    # assert (
    #     PROJECT_NAME != ""
    # ), "Project name is not set in the config file. Please set a project name"
    # ABRICATE = user_metadata["abricate"]
    #
    # PRIMER_FASTA = (
    #     user_metadata.get("primers_fasta")
    #     if user_metadata.get("primers_fasta") is not None
    #     else False
    # )
    # PRIMER_FASTA = (
    #     f"../user/primers/{user_metadata['primers_fasta']}" if PRIMER_FASTA else False
    # )
    #
    #
    # CONSENSUS_TOOL = user_metadata["illumina_consensus"]
    #
    # assert CONSENSUS_TOOL in "iVar snippy", "Consensus tool not set correctly"
    #
    # primers_setup = (
    #     PRIMER_FASTA is False or CONSENSUS_TOOL == "iVar" and os.path.isfile(PRIMER_FASTA)
    # )
    # if CONSENSUS_TOOL == "iVar":
    #     assert (
    #         primers_setup
    #     ), "Primers file not found. Please check the name and path to the primers file"
    #
    #
    # identification, version = get_identification_version(SEGMENTS, REFERENCE_GB)
    #
    # config_user = {
    #     "samples": samples,
    #     "project": user_metadata["project_name"],
    #     "locus": get_locus(REFERENCE_GB),
    #     "proteins": get_genes(REFERENCE_GB),
    #     "identification": identification,
    #     "version": version,
    #     "sample_type": sample_data.get_sample_type(),
    # }
    #
    # write_yaml(dic_directory["config_file"], config_user)
    # software_parameters = read_yaml(dic_directory["software_parameters"])
    #
    #
    #
