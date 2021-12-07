"""File for generating backbone swaps."""

from sgRNA_Generator import sgRNAGenerator
from Bio import SeqIO
import os
from Bio import pairwise2 as algi

alignments = algi.align.globalxx("AGTAGACCCC","AGAGACCATTGGGG")
print(algi.format_alignment(*alignments[0]))




class BackBoneSwap:

    def __init__(self):
        os.chdir("C:/Users/bmendoza/Documents/Auto_Ref_Test/bbswp/")

        # Parameters functions
        self.params_dict = dict()
        paramsfile = os.getcwd() + "/params.csv"
        self.params_parser(paramsfile)

        # Get the homarms:
        self.concatenated_arms = str()
        f = open(os.getcwd() + "/arms.csv")
        for line in f:
            myline = line[:-1].split(",")
            self.concatenated_arms += myline[1]
        f.close()

        # Containers for parts generated
        self.guides_and_arms = dict()

        for cluster in self.params_dict:
            self.get_guides(cluster)
            self.bbswap_primers_generator(cluster)
            for item in self.guides_and_arms[cluster]:
                print(cluster, item[0], item[1])

    # Imports the parameters from the parameter file
    def params_parser(self, pfile):
        f = open(pfile)
        for line in f:
            if not line.startswith("Cluster"):  # skip column title line
                newline = line[:-1].split("\t")
                self.params_dict[newline[0]] = newline[1:]
        f.close()

    def get_guides(self, cluster):
        sg = sgRNAGenerator("BBSWP", os.getcwd())
        self.guides_and_arms[cluster] = sg.guide_selector(cluster, self.params_dict[cluster][1:3])

    def bbswap_primers_generator(self, cluster):
        nt_seq = SeqIO.read(cluster + ".gb", "genbank").seq
        guide = self.guides_and_arms[cluster][0]
        homrev = nt_seq[guide[1]:guide[1]+30].reverse_complement()
        guide = self.guides_and_arms[cluster][1]
        homfwd = nt_seq[guide[1]-30:guide[1]]
        self.guides_and_arms[cluster].append((homrev, "rev"))
        self.guides_and_arms[cluster].append((homfwd, "fwd"))






