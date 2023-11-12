from scripts.yaml_io import read_yaml

class Abricate_Pangolin:
    def __init__(self, abricate_output, output, project):
        self.abricate_output = abricate_output
        self.output = output
        self.project = project

    def __call__(self, w):
        global checkpoints
        checkpoints.abricate_pangolin.get(**w)
        pattern = self.get_output()
        return pattern

    def get_output(self):
        abricate_dic = read_yaml(self.abricate_output)
        if abricate_dic["Species"] is None:
            abricate_dic["Species"] = []
        if abricate_dic["Genus"] is None:
            abricate_dic["Genus"] = []
        go_pangolin = bool(
            "coronavirus" in abricate_dic.get("Species")
            or "coronavirus" in abricate_dic.get("Genus")
        )
        pattern = (
            self.output
            if go_pangolin
            else expand(
                "projects/{project}/main_result/not_pangolin.csv", project=self.project
            )
        )
        return pattern
