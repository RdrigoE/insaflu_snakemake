import pandas
import yaml

class Data:
    def __init__(self,file):
        self.df = pandas.read_csv(file,na_values=None)
    
    def get_sample_names(self):
        return list(self.df["sample_name"])
    
    def get_sample_1(self):
        return list(self.df["fastq1"])
    
    def get_sample_2(self):
        return list(self.df["fastq2"])

    def get_sample_type(self):
        names = self.get_sample_names()
        type = self.df['tech']
        type_dic = {}
        for index in range(len(names)):
            type_dic[names[index]] = 'medaka' if type[index] == 'ont' else 'snippy'
        return type_dic


    def get_dic(self):
        dic = {}
        names = self.get_sample_names()
        for name in names:
            dic[name] = {}
        
        for idx,s1 in enumerate(self.get_sample_1()):
            
            dic[names[idx]]['fastq1'] = s1 if type(s1) == type('string') else None
        
        for idx,s2 in enumerate(self.get_sample_2()):
            dic[names[idx]]['fastq2'] = s2 if type(s2) == type('string') else None

        for idx,tech in enumerate(list(self.df['tech'])):
            dic[names[idx]]['tech'] = tech if type(tech) == type('string') else None
        
        return dic

    def get_options(self):
        project_dic = self.get_dic()
        final_dictionary = {}
        single_illumina={}
        paired_illumina={}
        ont_samples={}
        for names in project_dic:
            if project_dic[names]['fastq1'] and  project_dic[names]['fastq2'] and project_dic[names]['tech'] == 'illumina':
                paired_illumina[names] = project_dic[names]
                final_dictionary[names] = project_dic[names]
            elif project_dic[names]['fastq1'] and  not project_dic[names]['fastq2'] and project_dic[names]['tech'] == 'illumina':
                single_illumina[names] = project_dic[names]
                final_dictionary[names] = project_dic[names]
            elif project_dic[names]['fastq1'] and not project_dic[names]['fastq2'] and project_dic[names]['tech'] == 'ont':
                ont_samples[names] = project_dic[names]
                final_dictionary[names] = project_dic[names]
        return [paired_illumina,single_illumina,ont_samples, final_dictionary]
    def is_single_end(self):
        if self.get_sample_2() == []:
            return True
        return False
    
    def is_unique_file(self):
        if len(self.get_sample_1()) == 1 and self.get_sample_2() == []: 
            return True
        return False

def read_yaml(file):
    with open(file) as file:
        return yaml.load(file, Loader=yaml.FullLoader)


