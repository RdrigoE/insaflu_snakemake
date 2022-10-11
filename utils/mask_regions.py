from Bio import SeqIO
import Bio.Seq
from Bio.Seq import MutableSeq
import sys

def compute_masking_sites(sequence, ranges=None, singles_positions=None, from_beggining=None, from_end=None):
    length = len(sequence)
    masking_sites = []
    if ranges != None:
        for pos in ranges:
            pos[0] = int(pos[0])
            pos[1] = int(pos[1])
            if pos[1] < length and pos[1] > pos[0]:
                for position in range(pos[0]-1,pos[1]-1):
                    masking_sites.append(position)
    if single_positions != None:
        for pos in single_positions:
            pos = int(pos)
            if pos < length:
                masking_sites.append(pos-1)
    if from_beggining != None:
        for pos in range(int(from_beggining)):
            masking_sites.append(pos)
    if from_end != None:
        for pos in range(length - int(from_end), length):
            masking_sites.append(pos)
    return masking_sites

# consensus = '/home/reusebio/tese/insaflu_snakemake/projects/flu_testing_BVictoria/sample_SRR10885406/snippy/SRR10885406_consensus.fasta'
# single_positions = '1,2,10,400'.split(',')
# ranges = [x.split('-') for x in '10-30,40-50,60-90'.split(',')]
# from_beggining = '2'
# from_end = '2'

consensus = sys.argv[1]
output = sys.argv[2]

single_positions = sys.argv[3].split(',')
ranges = [x.split('-') for x in sys.argv[4].split(',')]
from_beggining = sys.argv[5]
from_end = sys.argv[6]

final_mask_consensus = []

for record in SeqIO.parse(consensus, "fasta"):
    sequence = MutableSeq(record.seq)
    masking_sites = compute_masking_sites(sequence, ranges, single_positions, from_beggining, from_end)
    ### Taken from insaflu
    ref_pos = 0
    ref_insertions = 0
    for _ in range(len(sequence)):
        if (sequence[_] == '-'):
            ref_insertions += 1
            continue
        if ref_pos in masking_sites:
            sequence[ref_pos + ref_insertions] = 'N'
        ref_pos += 1
        if ((ref_pos + ref_insertions) >= len(record.seq)): break
    ### End of insaflu code
    record.seq = sequence
    print(record.id, ' ===> ', record.seq)
    final_mask_consensus.append(record)

with open(output, "w") as handle_fasta_out:
	SeqIO.write(final_mask_consensus, handle_fasta_out, "fasta")