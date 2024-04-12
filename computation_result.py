import numpy as np

numpy_image_dict = {'gene=ORF1ab': [(115, 115), 7],
                    'gene=S': [(62, 62), 22],
                    'gene=ORF3a': [(28, 30), 12],
                    'gene=E': [(15, 16), 12], 
                    'gene=M': [(26, 27), 33],
                    'gene=ORF6': [(14, 14), 10],    
                    'gene=ORF7a': [(19, 20), 14],
                    'gene=ORF7b': [(12, 12), 12],
                    'gene=ORF8': [(19, 20), 14],
                    'gene=N': [(36, 36), 36],
                    'gene=ORF10': [(11, 11), 4]}

class DNA():
    def __init__(self, sequence):
        sequence = sequence.upper()
        sequence = sequence.replace(' ', '') 
        self.dir_3_5 = sequence
        self.dir_5_3 = self.dir_5_3_strand()
        self.mRNA = None
        self.num_array = None

    def transcription(self):
        trans = ''
        for nuc in self.dir_5_3:
            if nuc == 'A':
                trans += 'U'
            if nuc == 'T':
                trans += 'A'
            if nuc == 'C':
                trans += 'G'
            if nuc == 'G':
                trans += 'C'
            if nuc == 'N':
                trans += 'N'
        self.mRNA = trans
        return self.mRNA

    def dir_5_3_strand(self):
        dir_5_3 = ''
        for nuc in self.dir_3_5:
            if nuc == 'A':
                dir_5_3 += 'T'
            if nuc == 'T':
                dir_5_3 += 'A'
            if nuc == 'C':
                dir_5_3 += 'G'
            if nuc == 'G':
                dir_5_3 += 'C'
            if nuc == 'N':
                dir_5_3 += 'N'
        return dir_5_3

    def numpify(self):
        arr = ''
        for i in self.dir_3_5:
            if i == 'A':
                arr += '0 '
            if i == 'T':
                arr += '255 '
            if i == 'C':
                arr += '100 '
            if i == 'G':
                arr += '200 '
            if i == 'N':
                arr += '75 '   
        arr_np = np.fromstring(arr, dtype=np.uint8, sep=' ')        
        self.num_array = arr_np
        return self.num_array
    
def read_dna_seq(file_name):
    try:
        with open(file_name, 'r') as fil:
            fil_list = fil.readlines()
    except FileNotFoundError:
        try:
            with open('templates/sample1.txt', 'r') as fil:
                fil_list = fil.readlines()
        except FileNotFoundError:
            fil = open('reference1.txt', 'r')
            fil_list = fil.readlines()
    # fil = open(file_name, 'r')
    # fil_list = fil.readlines()
    fil.close
    
    genome = {}
    gene_name = ''
    gene_seq = ''
    for i in fil_list:
        if i[0] == '>':
            if list(genome.keys()) != []:
                gene_seq = gene_seq.replace('\n', '')
                genome[gene_name].append(gene_seq)
            gene_seq = ''
            g_st = i.find('[gene=')
            g_end = i[g_st:].find(']')

            if g_st > 0 and g_end > 0:
                gene_name = i[g_st+1:g_st+g_end]
                genome[gene_name] = []
        else:
            gene_seq += i
    gene_seq = gene_seq.replace('\n', '')
    genome[gene_name].append(gene_seq)
    return genome

def gene_mod(genome):
    genome_keys = list(genome.keys())
    for k in genome_keys:
        if len(numpy_image_dict[k]) > 1:
            N = numpy_image_dict[k][1]
            seq = add_N(N, genome[k][0])
            genome[k][0] = seq
    return genome

def add_N(n, seq):
    for i in range(0, n):
        seq += 'N'
    return seq

dict_seq_1 = read_dna_seq('sample1.txt')
dict_seq_1 = gene_mod(dict_seq_1)

dict_seq_2 = read_dna_seq('reference1.txt')
dict_seq_2 = gene_mod(dict_seq_2)

gene_names = list(numpy_image_dict.keys())
row = 0
col = 0
mut_dict = {}

for gene_name in gene_names:
    G = gene_name[5:]
    gene_a = DNA(dict_seq_1['gene=' + G][0])
    gene_a.transcription()
    numpify_a = gene_a.numpify()
    numpify_a = numpify_a.reshape(numpy_image_dict['gene=' + G][0])
    col += 1
    gene_b = DNA(dict_seq_2['gene=' + G][0])
    gene_b.transcription()
    numpify_b = gene_b.numpify()
    numpify_b = numpify_b.reshape(numpy_image_dict['gene=' + G][0])
    col += 1
    mut = numpify_b - numpify_a
    if mut.any():
        mut_nec = np.nonzero(mut)
        x = mut_nec[0]
        y = mut_nec[1]
        l = 0
        mut_dict[G] = []
        for i in x:
            us_base = numpify_a[i][y[l]]
            ch_base = numpify_b[i][y[l]]
            mut_base = mut[i][y[l]]
            info_list = [ch_base, us_base, mut_base, (i, y[l])]
            mut_dict[G].append(info_list)
            print("Mutated DNA Base {} in reference and Base {} in sample at position {} For the Gene {}" .format(ch_base, us_base, (i, y[l]), G))
            l += 1
    row += 1
    col = 0