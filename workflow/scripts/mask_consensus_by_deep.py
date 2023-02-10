"""tbd"""
import sys
from Bio import SeqIO

if __name__ == "__main__":
    reference_fasta = sys.argv[1]
    consensus_file = sys.argv[2]
    output_file = sys.argv[3]
    segment_name = sys.argv[4]

    with open(reference_fasta, "r", encoding="utf-8") as handle_fasta:
        dt_consensus = SeqIO.to_dict(SeqIO.parse(consensus_file, "fasta"))
        for record in SeqIO.parse(handle_fasta, "fasta"):
            if record.id == segment_name:  ### make mask
                ### get sequences
                vect_out_fasta_to_align = []
                record_id = record.id
                record.id = ""
                record.id = record_id + "_ref"
                # delete descriptions
                vect_out_fasta_to_align.append(record)
                vect_out_fasta_to_align.append(dt_consensus[record_id])
                with open(output_file, "w", encoding="utf-8") as handle_fasta_out_align:
                    SeqIO.write(
                        vect_out_fasta_to_align, handle_fasta_out_align, "fasta"
                    )
