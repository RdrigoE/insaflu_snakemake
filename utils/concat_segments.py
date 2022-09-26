import sys
from Bio import SeqIO
from get_locus import get_locus

def write_fasta(dictionary, filename):
        import textwrap
        with open(filename, "w") as fasta:
            for key, value in dictionary.items():
                fasta.write(f">{key}\n")
                fasta.write("\n".join(textwrap.wrap(str(value), 70)))
                fasta.write("\n")

alignment = sys.argv[1]
reference_gb = sys.argv[2]
species = sys.argv[3]
output = sys.argv[4]
counter = 0
set_read = 0
locus = get_locus(reference_gb,species)
new_file = []
reference_done = False


if type(locus) == type([0]):
    loop = len(locus)
    for record in SeqIO.parse(alignment, "fasta"):
        if counter == loop:
            counter = 0
            set_read += 1 
            if set_read == 1:
                reference_done = True
        if counter < loop:
            try:
                #print(f'Adding {counter} || Record Id: {record.id}')
                new_file[set_read].seq +=  record.seq
                # print(f'Added sequence {counter}')
            except:
                new_file.append(record)
                if reference_done:
                    end = record.id.index('__')
                    new_file[set_read].id = record.id[:end]
                else:
                    new_file[set_read].id = sys.argv[3]
                    print(new_file[set_read].id)
                print(reference_done)
        counter+=1
    dic = {}
    for identifier in new_file:
        dic[identifier.id] = identifier.seq

    write_fasta(dic, output)
    #SeqIO.write(new_file, output, "fasta")
#print(new_file[0].seq)
else:
    SeqIO.write(list(SeqIO.parse(alignment, "fasta")), output, "fasta")