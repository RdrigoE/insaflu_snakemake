import pandas


class Data:
    def __init__(self,file):
        self.df = pandas.read_csv(file,na_values=None)
    
    def get_sample_names(self):
        return list(self.df["sample_name"])
    
    def get_sample_1(self):
        return list(self.df["sample_1"])
    
    def get_sample_2(self):
        if self.df["sample_2"].isnull().values.any():
            return []
        return list(self.df["sample_2"])
    
    def is_single_end(self):
        if self.get_sample_2() == []:
            return True
        return False
    
    def is_unique_file(self):
        if len(self.get_sample_1()) == 1 and self.get_sample_2() == []: 
            return True
        return False