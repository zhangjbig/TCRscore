import os
from typing import NoReturn
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from itertools import combinations
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter


# Args parser
parser = ArgumentParser(description="Specifying Input Parameters")
parser.add_argument("-d", default="./hlahed/grantham_matrix.txt", help="Distance matrix for all amino acids, default: database/grantham_matrix.txt")
parser.add_argument("-f", default="./hlahed/ABC_prot.fa", help="Amino acid sequences in fasta format, default: database/ABC_prot.fa.")
parser.add_argument("-i", required=True, help="Input file with tab-delimited with individual HLA typing.")
parser.add_argument("-p", action="store_true", help="Paired HED score.")
parser.add_argument("-o", required=True, help="Output file.")
args = parser.parse_args()


def read_fasta(infile: str) -> dict:
    infile = Path(infile)
    record = SeqIO.parse(infile, "fasta")    ##  read fasta file
    sequences = {seq.id: str(seq.seq) for seq in record}
    seq_len = [len(value) for value in sequences.values()]
    if len(set(seq_len)) != 1:
        raise ValueError(f"Sequences have different lengths: {seq_len}. Please ensure all sequences are of the same length.")
    return sequences

def read_aa(infile: str) -> dict:
    infile = Path(infile)
    df = pd.read_csv(infile, header=0, sep="\t", index_col=0)
    aa_pairwise_dis = df.to_dict()
    return aa_pairwise_dis

def calculate_distance(hla1: str, hla2: str, sequences: dict, distance: dict) -> float:
    seq_hla1 = sequences.get(hla1, None)
    seq_hla2 = sequences.get(hla2, None)
    if not seq_hla1 or not seq_hla2:
        return("NA")
    else:
        seq_len = len(seq_hla1)
        dis = 0
        for i in range(seq_len):
            aa1 = seq_hla1[i]
            aa2 = seq_hla2[i]
            try:
                dis += distance[aa1][aa2]
            except KeyError:
                raise ValueError(f"Amino acid pair ({aa1}, {aa2}) not found in the distance matrix.")
        dis = dis/seq_len
        return dis


def main():
    sequences = read_fasta(args.f)
    aa_pairwise_dis = read_aa(args.d)
    infile = Path(args.i)
    outfile = Path(args.o)

    df = pd.read_csv(infile, header=0, sep="\t")
    if args.p:
        df2 = pd.melt(df, id_vars=["Sample"], value_vars=["A1", "A2", "B1","B2", "C1","C2"])
        alleles = set(df2["value"].values.tolist())
        alleles_pair = combinations(alleles, 2)
    
        outheader = ["Allele1","Allele2","HED"]
        with open(outfile, "w") as fw:
            fw.write("\t".join(outheader) + "\n")
            for allele1, allele2 in alleles_pair:
                dis_hla_pair = calculate_distance(allele1, allele2, sequences, aa_pairwise_dis)
                outline = [allele1, allele2, dis_hla_pair]
                outline = [str(x) for x in outline]

                fw.write("\t".join(outline) + "\n")
    else:
        outheader = ["Sample","HED_A","HED_B","HED_C","Mean_HE"]
        with open(outfile, "w") as fw:
            fw.write("\t".join(outheader) + "\n")
            for _, line in df.iterrows():
                hla_a1 = line["A1"]
                hla_a2 = line["A2"]
                dis_hla_a = calculate_distance(hla_a1, hla_a2, sequences, aa_pairwise_dis)

                hla_b1 = line["B1"]
                hla_b2 = line["B2"]
                dis_hla_b = calculate_distance(hla_b1, hla_b2, sequences, aa_pairwise_dis)
                
                hla_c1 = line["C1"]
                hla_c2 = line["C2"]
                dis_hla_c = calculate_distance(hla_c1, hla_c2, sequences, aa_pairwise_dis)

                if dis_hla_a == "NA" or dis_hla_b == "NA" or dis_hla_c == "NA":
                    dis_mean = "NA"
                else:
                    dis_mean = (dis_hla_a + dis_hla_b + dis_hla_c) / 3

                outline = [line["Sample"], dis_hla_a, dis_hla_b, dis_hla_c, dis_mean]
                outline = [str(x) for x in outline]

                fw.write("\t".join(outline) + "\n")

if __name__ == "__main__":
    main()