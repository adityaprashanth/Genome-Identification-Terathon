from flask import *
import numpy as np

app = Flask(__name__) 

@app.route('/') 
def main(): 
    return render_template("home.html")

@app.route('/success', methods=['POST']) 
def success(): 
    f = request.files['file'] 
    f.save(f.filename) 
    if request.method == 'POST': 
        soln = ""

        numpy_image_dict = {'gene=1080':[(78,78),14],
                            'gene=3077':[(33,33),6],
                            'gene=3043':[(25,26),22],}

        class DNA:
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
            fil = open(file_name, 'r')
            fil_list = fil.readlines()
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

        dict_seq_1 = read_dna_seq('the_sample.txt')
        dict_seq_1 = gene_mod(dict_seq_1)

        dict_seq_2 = read_dna_seq('the_reference.txt')
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
                    soln += ("\nMutated DNA Base {} in reference and Base {} in sample at position {} For the Gene {}\n".format(ch_base, us_base, (i, y[l]), G))
                    if(G == "1080"):
                        soln += ("chances of catching cystic fibrosis, if not at least non symptomatic carrier\n");
                    if(G == "3077"):
                        soln += ("chances of catching hemochromatosis, if not at least non symptomatic carrier\n");
                    if(G == "3043"):
                        soln += ("chances of catching sickle cell anemia, if not at least nonn symptomatic carrier\n");
                    l += 1
            row += 1
            col = 0
        return render_template("Acknowledgement.html", name = soln) 

@app.route('/signup', methods=['GET', 'POST'])  # Add route for signup page
def signup():
    if request.method == 'POST':
        # Handle signup form submission
        return redirect('/')
    else:
        return render_template("signup.html")

if __name__ == '__main__':
    app.run(debug=True)
