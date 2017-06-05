import collections

def protein_translate(rna, codons_table):
    res = ""

    for i in range(0, len(rna), 3):
        if len(rna[i:]) < 3:
            break
        if codons_table[rna[i: i + 3]] is "*":
            break
        res += codons_table[rna[i: i + 3]]

    return res

def rna_encodes(rna, prot, codons_table):
    res = list()
    pot = 3 * len(prot)
    for i in range(0, len(rna) - pot + 1):
        if protein_translate(rna[i: i + 3*len(prot)], codons_table) == prot:
            res.append(rna[i:i+len(prot)*3])
        if protein_translate(reverse_comp_rna(rna[i: i + 3*len(prot)]), codons_table) == prot:
            res.append(rna[i:i+len(prot)*3])

    return res

##########################################################################

def parse_input_codons_table(fname):
    codons_table = dict()

    workfile = open(fname, 'r')
    codons = workfile.readlines()
    codons_pair = map(lambda x: x.strip(), codons)

    for pair in codons_pair:
        (codon, amin) = pair.split(" ")
        codons_table[codon] = amin

    return codons_table

def parse_input_rna(fname):
    workfile = open(fname, 'r')
    rna = workfile.readline().strip()
    return rna

def parse_input_dna_and_prot(fname):
    workfile = open(fname, 'r')
    rna = workfile.readline().strip()
    prot = workfile.readline().strip()
    return (rna, prot)

def dna_to_rna(dna):
    def t2u(c):
        if c is "T":
            return "U"
        else:
            return c

    return "".join(t2u(c) for c in dna)

def rna_to_dna(dna):
    def u2t(c):
        if c is "U":
            return "T"
        else:
            return c

    return "".join(u2t(c) for c in dna)

def reverse_comp_dna(dna):
    def reverse_char(c):
        if c is "A":
            return "T"
        elif c is "T":
            return "A"
        elif c is "G":
            return "C"
        else:
            return "G"
    res = "".join([reverse_char(c) for c in dna][::-1])
    return res

def reverse_comp_rna(dna):
    def reverse_char(c):
        if c is "A":
            return "U"
        elif c is "U":
            return "A"
        elif c is "G":
            return "C"
        else:
            return "G"
    res = "".join([reverse_char(c) for c in dna][::-1])
    return res

def format_list(l):
    for e in l:
        print(e)

##########################################################################

def main():
    codons = parse_input_codons_table("codons.txt")
    (dna, prot) = parse_input_dna_and_prot("dataset_96_7.txt")
    res = [rna_to_dna(x) for x in rna_encodes(dna_to_rna(dna), prot, codons)]
    format_list(res)

if __name__ == '__main__':
    main()
