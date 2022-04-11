import sys


compliment = {'A':'T', 'T':'A', 'C': 'G', 'G':'C'}
START_CODONS = ['ATG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']
translatinon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

def lost_in_translation(ans_nucl):
    ans=[]
    for orf in ans_nucl:
        ans1=[]
        for j in range(0, len(orf), 3):
            codon=orf[j:j+3]
            ans1.append(translatinon_table[codon])
        ans.append("".join(ans1))

    return ans

def find_ORF(seq):
    start_codons = []
    stop_codons = []
    ans_nucl = []
    start_codons=find_start_codon(seq, start_codons)
    stop_codons=find_stop_codon(seq, stop_codons)
    for i in start_codons:
        for j in stop_codons:
            if (i - j) % 3 == 0 and i<j:
                ans_nucl.append(seq[i:j + 3])
    return ans_nucl


def find_start_codon(seq, start_codons):

    for i in range(0, len(seq)-3, 1):
        if seq[i:i+3] in START_CODONS:
            start_codons.append(i)

    return start_codons


def find_stop_codon(seq, stop_codons):
    j = 0
    for i in range(0, len(seq)-3 , 1):
        if seq[i:i+3] in STOP_CODONS:
            stop_codons.append(i)
    return stop_codons

def make_reverse(seq):
    return seq[::-1]


def make_compliment(seq):
    seq=make_reverse(seq)
    return "".join([compliment[x] for x in seq])

def read_fasta(file):
    m=file.readlines()[1:]
    fasta_str = "".join(m)
    fasta_str=fasta_str.replace("\n" , "")
    fasta_str=fasta_str.replace("\r", "")
    return fasta_str


if __name__ == '__main__':
    fin = open('srt.txt' , 'r')
    seq=read_fasta(fin)
    seq_rev=make_compliment(seq)
    l=find_ORF(seq)
    m=find_ORF(seq_rev)
    acids=lost_in_translation(l)
    f = open('output.txt', 'w')
    for element in acids:
        f.write(element)
        f.write('\n')
    f.close()
    acid_2=lost_in_translation(m)



