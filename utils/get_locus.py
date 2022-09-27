from Bio import SeqIO

def get_locus(genbank_file,possible_name):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name. If there are multiple records, 
    it will return a list of all locus numbers.
    
    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    count = 1
    locus = []
    handle_gb = open(genbank_file)
    for record in SeqIO.parse(handle_gb, "genbank"):
        locus.append(count)
        count+=1
    handle_gb.close()
    if len(locus) == 1:
        return possible_name
    else:
        return locus
