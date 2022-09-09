# with open('config_user/config_run.yaml') as file:
#     config_user = yaml.load(file, Loader=yaml.FullLoader)

# rule pre_snpeff:
#     input:
#         ref = REFERENCE

#     output:
#         directory("projects/{project}/main_result/snpeff/genomes")   
#     conda:
#         "../envs/snpeff.yaml"
#     shell:
#         "mkdir -p {output} |"
#         "cp {input.ref} projects/{wildcards.project}/main_result/snpeff/genomes/"


# rule snpeff:
#     input:
#         i = "projects/{project}/main_result/freebayes/{sample}_var.vcf",
#         ref = REFERENCE_GFF3
#     output:
#         o = "projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",
#     conda:
#         "../envs/snpeff.yaml"
#     params:
#         "-no-downstream -no-upstream -no-intergenic -no-utr -noStats -c config/snpeff.config "
#     shell:
#         "snpEff {params} -v {input.ref} {input.i}  > {output.o}"


# def run_snpEff(self, fasta_file, genbank, vcf_file, out_file):
# 		"""
# 		./snpEff ann -no-downstream -no-upstream -no-intergenic -no-utr -c ../path_to/reference/snpeff.config -dataDir . -noStats ref sample.vcf > sample_annot.vcf
# 		"""
# 		temp_dir = self.utils.get_temp_dir()
# 		(fasta_file_name, snpeff_config) = self.get_snpeff_config(fasta_file)
		
# 		temp_vcf_file = os.path.join(temp_dir, os.path.basename(vcf_file))
# 		self.utils.copy_file(vcf_file, temp_vcf_file)
		
# 		## create the database
# 		out_gff_file = self.utils.get_temp_file('temp_gff', '.gff')
# 		self.run_genbank2gff3(genbank, out_gff_file)
		
# 		### count sequences, if none return None
# 		if (self.utils.get_number_sequeces_in_gff_file(out_gff_file) == 0):
# 			os.unlink(out_gff_file)
# 			os.unlink(snpeff_config)
# 			self.utils.remove_dir(temp_dir)
# 			return None

# 		datase_dir = "{}".format(fasta_file_name)
# 		cmd = "mkdir -p {}".format(os.path.join(temp_dir, datase_dir)); os.system(cmd)
# 		cmd = "mkdir -p {}".format(os.path.join(temp_dir, 'genomes')); os.system(cmd)
# 		self.utils.copy_file(out_gff_file, os.path.join(temp_dir, datase_dir, 'genes.gff'))
# 		temp_file = os.path.join(temp_dir, 'genomes', fasta_file_name + '.fa')
# 		self.utils.copy_file(fasta_file, temp_file)
# 		os.unlink(out_gff_file)
		
# 		## indexing database
# 		## snpEff build -c reference/snpeff.config -dataDir . -gff3 ref 2>> run_snippy2_1.log
# 		cmd = "%s build -c %s -dataDir %s -gff3 %s" % (self.software_names.get_snp_eff(),\
# 						snpeff_config, temp_dir, fasta_file_name)
# 		exist_status = os.system(cmd)
# 		if (exist_status != 0):
# 			os.unlink(snpeff_config)
# 			self.logger_production.error('Fail to run: ' + cmd)
# 			self.logger_debug.error('Fail to run: ' + cmd)
# 			raise Exception("Fail to create snpEff database")
		
# 		### create the annotation
# 		cmd = "%s ann %s -c %s -dataDir %s %s %s > %s" % (self.software_names.get_snp_eff(), self.software_names.get_snp_eff_parameters(),\
# 						snpeff_config, temp_dir, fasta_file_name, temp_vcf_file, out_file)
# 		exist_status = os.system(cmd)
# 		if (exist_status != 0):
# 			os.unlink(snpeff_config)
# 			self.logger_production.error('Fail to run: ' + cmd)
# 			self.logger_debug.error('Fail to run: ' + cmd)
# 			raise Exception("Fail to run snpEff")
		
# 		#### add the transform p.Val423Glu to p.V423G
# 		parse_out_files = ParseOutFiles()
# 		out_file_transformed_amino = parse_out_files.add_amino_single_letter_code(out_file)
# 		self.utils.move_file(out_file_transformed_amino, out_file)
		
# 		os.unlink(snpeff_config)
# 		self.utils.remove_dir(temp_dir)
# 		return out_file