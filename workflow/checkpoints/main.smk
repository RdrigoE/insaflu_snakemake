import csv

from scripts.extract_gb_info import get_locus_and_genes


class Checkpoint_Main:
    def __init__(
        self,
        output_tuple,
        genbank_file,
        coverage_file,
        coverage_limit,
        project,
        samples,
        segments
    ):
        self.output_tuple = output_tuple
        self.genbank_file = genbank_file
        self.coverage_file = coverage_file
        self.coverage_limit = coverage_limit
        self.project = project
        self.samples = samples
        self.segments = segments
    def __call__(self, w):
        global checkpoints

        checkpoints.mergeCoverage.get(**w)

        pattern = self.get_output()
        pattern.extend(self.get_output_prefix_sufix())

        if pattern == []:
            return "template_empty.txt"
        return pattern

    def get_locus_w_coverage_prefix_sufix(self):
        with open(self.coverage_file, newline="") as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=",")
            coverage_list = list(csv_reader)
            coverage_list = list(map(lambda x: x[1:], coverage_list))
            for idx, li in enumerate(coverage_list):
                coverage_list[idx] = list(map(lambda x: float(x), li))

            final_output = []
            for sample in coverage_list:
                if min(sample) >= self.coverage_limit:
                    final_output.append(True)
                else:
                    final_output.append(False)

        return final_output

    def get_output_prefix_sufix(self):
        locus_protein = []
        valide_locus = self.get_locus_w_coverage_prefix_sufix()
        segments = get_locus_and_genes(self.genbank_file)
        if True in valide_locus:
            for prefix, sufix in self.output_tuple:
                for idx, seg in enumerate(segments, start=0):
                    for gene in segments[seg]:
                        locus_protein.append(
                            f"{prefix}{seg}/Alignment_aa_{seg}_{gene}{sufix}"
                        )
        return locus_protein

    def get_locus_w_coverage(self):
        with open(self.coverage_file, newline="") as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=",")
            coverage_list = list(csv_reader)
            coverage_list = list(map(lambda x: x[1:], coverage_list))
            for idx, li in enumerate(coverage_list):
                coverage_list[idx] = list(map(lambda x: float(x), li))

            final_output = []
            for sample in coverage_list:
                if min(sample) >= self.coverage_limit:
                    final_output.append(True)
                else:
                    final_output.append(False)

        return final_output

    def get_output(self):
        return_list = []
        valide_locus = self.get_locus_w_coverage()

        if True in valide_locus:
            files = [
                expand(
                    "projects/{project}/main_result/validated_variants.csv",
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/Alignment_nt_All.fasta",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/All_nt_only_90plus.fasta",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/AllConsensus.fasta",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/All_nt.fasta",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/All_nt.nex",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/AllConsensus.nex",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/Alignment_nt_All.nex",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/All_nt_only_90plus.nex",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/snp_ready.txt",
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/Tree_ML_All.tree",
                    sample=self.samples,
                    project=self.project,
                ),
                expand(
                    "projects/{project}/main_result/Tree_ML_{seg}.tree",
                    sample=self.samples,
                    project=self.project,
                    seg=self.segments,
                ),
            ]
            if user_metadata["get_minor_variants"]:
                files.extend(
                    [
                        expand(
                            "projects/{project}/main_result/validated_minor_iSNVs.csv",
                            project=self.project,
                        ),
                        expand(
                            "projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",
                            project=self.project,
                        ),
                        expand(
                            "projects/{project}/main_result/proportions_iSNVs_graph.csv",
                            project=self.project,
                        ),
                    ]
                )
            return_list = []
            for new in files:
                if isinstance(new, list):
                    for new_2 in new:
                        return_list.append(new_2)
                else:
                    return_list.append(new)
            return return_list
        else:
            return expand("projects/{p}/main_result/warning.txt", p=self.project)
