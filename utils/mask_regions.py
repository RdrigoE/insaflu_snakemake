from Bio import SeqIO
import Bio.Seq
from Bio.Seq import MutableSeq
import sys
import argparse

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

parser = argparse.ArgumentParser()
parser.add_argument("consensus", type=str,
                    help="Consensus File Path")
parser.add_argument("output", type=str,
                    help="Output File Path")
parser.add_argument("-r","--regions", type=str,
                    help="Pass a string with the format x-y,n-z,...")
parser.add_argument("-s", "--single", type=str,
                    help="Pass a string with the format x,y,n,z,...")
parser.add_argument("-b", "--beggining", type=int,
                    help="Pass a integer with the format x")
parser.add_argument("-e", "--end", type=int,
                    help="Pass a integer with the format x")

args = parser.parse_args()

consensus = None
output = None
single_positions = None
ranges = None
from_beggining = None
from_end = None


if args.consensus: 
    consensus = args.consensus
if args.output: 
    output = args.output
if args.single: 
    single_positions = args.single.split(',')
if args.regions: 
    ranges = [x.split('-') for x in args.regions.split(',')]
if args.beggining: 
    from_beggining = args.beggining
if args.end: 
    from_end = args.end

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
    final_mask_consensus.append(record)

with open(output, "w") as handle_fasta_out:
	SeqIO.write(final_mask_consensus, handle_fasta_out, "fasta")