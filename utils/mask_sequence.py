def mask_sequence(self, sequence_ref, sequence_consensus, mask_sites,
					mask_from_beginning, mask_from_end, mask_range):
		"""Mask characters at the given sites in a single sequence record, modifying the
		record in place.
		Parameters
		----------
		sequence_ref : Bio.SeqIO.SeqRecord
		    A sequence that act like a reference. All the positions are referent to this seq.
		sequence_consensus : Bio.SeqIO.SeqRecord
		    A sequence to be masked
		mask_sites: list[int]
		    A list of site indexes to exclude from the FASTA.
		mask_from_beginning: int
		    Number of sites to mask from the beginning of each sequence (default 0)
		mask_from_end: int
		    Number of sites to mask from the end of each sequence (default 0)
		mask_invalid: bool
		    Mask invalid nucleotides (default False)
		Returns
		-------
		Bio.SeqIO.SeqRecord
		    Masked sequence in its original record object
		"""
		if len(sequence_ref.seq) == 0 or len(sequence_consensus.seq) == 0: return sequence_consensus
		# need to align
		sequence_length = len(sequence_ref.seq)
		seq_ref, seq_consensus = self.align_two_sequences(str(sequence_ref.seq), str(sequence_consensus.seq))
		
		# Convert to a mutable sequence to enable masking with Ns.
		beginning = int(mask_from_beginning) if not mask_from_beginning is None and len(mask_from_beginning) > 0 else -1 
		end = int(mask_from_end) if not mask_from_end is None and len(mask_from_end) > 0 else -1
		
		## collecting all positions to maks
		dt_positions = {}
		if beginning != -1:
			for _ in range(0, beginning): dt_positions[_] = 1
		if end != -1:
			if ( (len(str(sequence_ref.seq)) - end) < 0): pos_from = 0
			else: pos_from = len(str(sequence_ref.seq)) - end
			for _ in range(pos_from, len(str(sequence_ref.seq))): dt_positions[_] = 1

		## several sites
		if not mask_sites is None and  len(mask_sites.split(',')[0]) > 0:
			for site in [int(_) - 1 for _ in mask_sites.split(',')]:
				if site < sequence_length: dt_positions[site] = 1
		## several ranges
		if not mask_range is None:
			for data_ in mask_range.split(','):
				if (len(data_) > 0 and len(data_.split('-')) == 2):
					for site in range(int(data_.split('-')[0]) - 1, int(data_.split('-')[1])):
						if site < sequence_length: dt_positions[site] = 1
		
		## mask positions
		masked_sequence = MutableSeq(seq_consensus)
		ref_insertions = 0
		ref_pos = 0
		for _ in range(len(seq_ref)):
			if (seq_ref[_] == '-'):
				ref_insertions += 1
				continue
			if ref_pos in dt_positions:
				masked_sequence[ref_pos + ref_insertions] = 'N'
			ref_pos += 1
			if ((ref_pos + ref_insertions) >= len(seq_consensus)): break
			
		sequence_consensus.seq = Seq(str(masked_sequence).replace('-', ''))
		return sequence_consensus
	
	def mask_sequence_by_sites(self, reference_fasta_file, consensus_fasta_file, genetic_elemets):
		""" masking consensus file with positions related with reference elements  """
		vect_record_out = []
		## always work with the backup	
		with open(reference_fasta_file, "rU") as handle_ref, open(consensus_fasta_file, "rU") as handle_consensus:
			dict_record_ref = SeqIO.to_dict(SeqIO.parse(handle_ref, "fasta"))
			for record_consensus in SeqIO.parse(handle_consensus, "fasta"):
				masking_consensus = genetic_elemets.dt_elements_mask.get(record_consensus.id, MaskingConsensus())
				if masking_consensus.has_data() and record_consensus.id in dict_record_ref:
					vect_record_out.append(self.mask_sequence(
						dict_record_ref[record_consensus.id], record_consensus,
						masking_consensus.mask_sites, masking_consensus.mask_from_beginning, 
						masking_consensus.mask_from_ends, masking_consensus.mask_regions))
				else: vect_record_out.append(record_consensus)

		if (len(vect_record_out) > 0):
			temp_file = self.utils.get_temp_file("masked_seq_", ".fasta")
			with open(temp_file, "w") as handle_fasta_out:
				SeqIO.write(vect_record_out, handle_fasta_out, "fasta")

			### move temp consensus to original position, if has info
			if os.stat(temp_file).st_size > 0:
				self.utils.move_file(temp_file, consensus_fasta_file)
			else: os.unlink(temp_file)