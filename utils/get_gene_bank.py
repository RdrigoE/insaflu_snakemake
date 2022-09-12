from Bio import SeqIO

genbank_file = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb"

def get_positions_gb(genbank_file):
    positions = []

    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        for features in record.features:
            if (features.type == 'CDS'):
                positions.append([features.qualifiers["gene"][0], features.location])

    positions_clean = []

    for idx,gene in enumerate(positions):
        positions_clean.append([gene[0],[]])
        for part in gene[1].parts:
            positions_clean[idx][1].append([int(part.start),int(part.end)]) 
    handle_gb.close()
    return positions_clean


def get_genes(genbank_file):
    genes = []

    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        for features in record.features:
            if (features.type == 'CDS'):
                genes.append(features.qualifiers["gene"][0])
    handle_gb.close()

    return genes

