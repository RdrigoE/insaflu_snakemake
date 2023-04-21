from workflow.scripts.add_freq_medaka_consensus import get_type_variation


class Test_get_type_variantion:
    def test_snp(self):
        assert get_type_variation(1, 1) == "snp"

    def test_insertion(self):
        assert get_type_variation(2, 3) == "ins"

    def test_deletion(self):
        assert get_type_variation(2, 1) == "del"

    def test_mnp(self):
        assert get_type_variation(2, 2) == "mnp"

    def test_complex(self):
        assert get_type_variation(0, 0) == "complex"
