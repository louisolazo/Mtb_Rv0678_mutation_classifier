import argparse
import re
import random

# Global variables
REF_GENE = "AGCGTCAACGACGGGGTCGATCAGATGGGCGCCGAGCCCGACATCATGGAATTCGTCGAACAGATGGGCGGCTATTTCGAGTCCAGGAGTTTGACTCGGTTGGCGGGTCGATTGTTGGGCTGGCTGCTGGTGTGTGATCCCGAGCGGCAGTCCTCGGAGGAACTGGCGACGGCGCTGGCGGCCAGCAGCGGGGGGATCAGCACCAATGCCCGGATGCTGATCCAATTTGGGTTCATTGAGCGGCTCGCGGTCGCCGGGGATCGGCGCACCTATTTCCGGTTGCGGCCCAACGCTTTCGCGGCTGGCGAGCGTGAACGCATCCGGGCAATGGCCGAACTGCAGGACCTGGCTGACGTGGGGCTGAGGGCGCTGGGCGACGCCCCGCCGCAGCGAAGCCGACGGCTGCGGGAGATGCGGGATCTGTTGGCATATATGGAGAACGTCGTCTCCGACGCCCTGGGGCGATACAGCCAGCGAACCGGAGAGGACGACTGA"
GENE_START = 778990
GENE_END = 779487
BP_LENGTH = len(REF_GENE)


def main():
    parser = argparse.ArgumentParser(description="MTB Rv0678 Mutation Classifier")
    parser.add_argument("--snp", type=str, help='Input SNP string in standard notation (eg. "c.492C>T" or "g.779440G>T"). Accepts CDS (c.) and genomic (g.) notations. Wrap notation with double-quotes(i.e. "<notation>").', required=False)
    args = parser.parse_args()

    snp = args.snp

    if not snp:
        print(f"\nGenerating random SNP...")
        random_pos = int(random.randrange(BP_LENGTH))
        base_options = [i for i in "ACGT" if i != (REF_GENE[random_pos])]
        mut_base = random.choice(base_options)
        snp = f"c.{random_pos+1}{REF_GENE[random_pos]}>{mut_base}"

    match = re.match(r"^(c|g)\.(\d+)([ACGT])>([ACGT])$", snp)
    if not match:
        print("Invalid SNP format.")
        return

    notation, pos, ref, mut = match.groups()
    pos = int(pos)

    if notation == "c":
        codon_no, ref_codon, mut_codon = cds_parser(pos, ref, mut)
    else:
        codon_no, ref_codon, mut_codon = genomic_parser(pos, ref, mut)



    ref_aa = translate_codon(ref_codon)
    mut_aa = translate_codon(mut_codon)
    mutation_type = classify_mutation(ref_aa, mut_aa)

    print (f"\nMycobacterium tuberculosis Rv0678 Mutation Classifier v0.1 \n")
    print (f"Input SNP: {snp}")
    print (f"Codon No.: {codon_no}")
    print (f"Codon: {ref_codon} -> {mut_codon}")
    print (f"Amino Acid: {ref_aa} -> {mut_aa}")
    print (f"Mutation Type: {mutation_type}\n")


def cds_parser(pos, ref, mut):

    coordinate = pos - 1
    if coordinate < 0 or coordinate >= BP_LENGTH:
        raise ValueError("Position out of CDS range.")

    if ref != REF_GENE[coordinate]:
        raise ValueError(f"Reference base {ref} does not match at position {pos}.")

    if mut not in "ACGT":
        raise ValueError(f"Invalid mutation base {mut}. Must be A, C, G, or T.")

    codon_start = (coordinate // 3) * 3
    ref_codon = REF_GENE[codon_start:codon_start + 3]

    if len(ref_codon) != 3:
        raise ValueError("Position does not align with a complete codon.")

    codon_pos = coordinate % 3
    mut_codon = ref_codon[:codon_pos] + mut + ref_codon[codon_pos + 1:]

    codon_no = coordinate // 3 + 1

    return codon_no, ref_codon, mut_codon


def genomic_parser(pos, ref, mut):

    if pos < GENE_START or pos > GENE_END:
        raise ValueError("Position out of gene range.")

    coordinate = pos - GENE_START

    if ref != REF_GENE[coordinate]:
        raise ValueError(f"Reference base {ref} does not match at position {pos}.")

    if mut not in "ACGT":
        raise ValueError(f"Invalid mutation base {mut}. Must be A, C, G, or T.")


    codon_start = (coordinate // 3) * 3
    ref_codon = REF_GENE[codon_start:codon_start + 3]

    if len(ref_codon) != 3:
        raise ValueError("Position does not align with a complete codon.")

    codon_pos = coordinate % 3
    mut_codon = ref_codon[:codon_pos] + mut + ref_codon[codon_pos + 1:]

    codon_no = coordinate // 3 + 1

    return codon_no, ref_codon, mut_codon



def translate_codon(codon):

    codon_table = {
    'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',    # Serine
    'TTC': 'Phe', 'TTT': 'Phe',    # Phenylalanine
    'TTA': 'Leu', 'TTG': 'Leu',    # Leucine
    'TAC': 'Tyr', 'TAT': 'Tyr',    # Tirosine
    'TAA': '*', 'TAG': '*',    # Stop
    'TGC': 'Cys', 'TGT': 'Cys',    # Cisteine
    'TGA': '*',    # Stop
    'TGG': 'Trp',    # Tryptofan
    'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',    # Leucine
    'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',    # Proline
    'CAC': 'His', 'CAT': 'His',    # Histidine
    'CAA': 'Gln', 'CAG': 'Gln',    # Glutamine
    'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',    # Arginine
    'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile',    # Isoleucine
    'ATG': 'Met',    # Methionine
    'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',    # Threonine
    'AAC': 'Asn', 'AAT': 'Asn',    # Asparagine
    'AAA': 'Lys', 'AAG': 'Lys',    # Lysine
    'AGC': 'Ser', 'AGT': 'Ser',    # Serine
    'AGA': 'Arg', 'AGG': 'Arg',    # Arginine
    'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',    # Valine
    'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',    # Alanine
    'GAC': 'Asp', 'GAT': 'Asp',    # Aspartic Acid
    'GAA': 'Glu', 'GAG': 'Glu',    # Glutamic Acid
    'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly'     # Glycine
}

    return codon_table.get(codon)


def classify_mutation(ref_aa, mut_aa):
    if ref_aa == mut_aa:
        return "Silent"
    elif mut_aa == "?":
        return "Unknown"
    elif ref_aa == "*":
        return "Stop Gain"
    elif mut_aa == "*":
        return "Stop Loss"
    else:
        return "Missense"



if __name__ == "__main__":
    main()
