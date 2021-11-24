"""File for generating PROSWP primers and sgRNA templates"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
from sgRNA_Generator import sgRNAGenerator

class AutoRefactor:

    def __init__(self, workingdir, mod_style, marker=True):
        os.chdir(workingdir)
        self.mod_style = mod_style

        """Dictionary containing relevant information for clusters including:
            -Name of Cluster
            -Paradigm of Refactoring (Resistance, Core, Gene Maximal)
            -CDS's of the cluster in order of priority to refactor
            -CDS's of clusters that need to be deleted (flagged as DEL- in name)"""
        self.params_dict = dict()
        self.promoters = dict()  # list of promoters that can be used with accompanying strengths in their tuples

        # Function to fill the above objects
        self.import_csv_marker(workingdir)

        # Initiate the sgRNA generator class object for getting guide sequences
        self.sgen = sgRNAGenerator("PROSWP", workingdir)

        # Now go through the clusters that have been imported:
        output_dict = dict()
        for cluster in self.params_dict:
            output_dict[cluster] = list()
            for modification in self.params_dict[cluster]:
                output_dict[cluster].append(self.sgen.guide_selector(cluster, modification))
            output_dict[cluster].append(self.hifi_hom_arm_generator(cluster))
            self.output(output_dict)

    """Function opens the output stream to the file for building the factory order."""
    def output(self, out):
        out_file_name = "myoutput.csv"
        f = open(out_file_name, 'w')
        for cluster in out:
            f.write(cluster)
            # This separates the refactoring site list
            for site in out[cluster]:
                # this separates the guides from the hom arms
                for item in site:
                    newline = str()
                    # and this separates the sequences and locations
                    for thing in item:
                        newline += ("," + str(thing))
                    f.write(newline + "\n")
        f.close()

    """Import the csv file of parameters into a csv that defines which operons to refactor."""
    def import_csv_marker(self, filedir):
        os.chdir(filedir)
        params_file_name = "proswp/params.csv"
        f = open(params_file_name)
        skipfirstline = True

        for line in f:
            if skipfirstline:
                skipfirstline = False
                continue
            params = line[:-1].split(",")
            # If new fasta_id generate a new dictionary item
            if params[0] in self.params_dict:
                self.params_dict[params[0]].append(params[1:])
            # Otherwise create a new list of parameters
            else:
                self.params_dict[params[0]] = [params[1:]]
        f.close()

        # And to import the promoter sequences:
        # Note: download sequences from LIMS and promoter strength as attribute later
        y = open("LIMS_DNA_parts.csv")
        for line in y:
            prom = line[:-1].split(",")
            self.promoters[prom[0]] = prom[1]
        y.close()

    def import_clusters_fasta(self):
        os.chdir(os.curdir + "/Refactoring_Clusters/")
        for file in os.listdir():
            if file.endswith(".fa") or file.endswith(".fasta"):
                cluster_seq = ""
                cluster_name = ""
                f = open(file)
                for line in f:
                    if line.startswith(">"):
                        cluster_name = line[:-1]
                    else:
                        cluster_seq += line[:-1]
                f.close()

    # Placeholder function for parsing .gb files
    def parse_clusters(self,filedir):
        os.chdir(filedir)
        for cluster in self.params_dict:
            id_dict = SeqIO.to_dict(SeqIO.parse(cluster + ".gb", "gb"))
            #for CDS_to_refactor in self.params_dict[cluster]:
                #id_dict[CDS_to_refactor].location

    """Inputs: The parameters and parts files, the genbank of the cluster.
        Output: The primer sequences for generating the homology arms."""
    def hifi_hom_arm_generator(self, cluster):
        output_list = list()
        nt_seq = SeqIO.read(cluster + ".gb", "genbank").seq
        for modification in self.params_dict[cluster]:
            # Get the overlap regions on the cluster
            print(modification)
            hom_fwd = nt_seq[int(modification[3])-30:int(modification[3])]  # hardcoded index needs to be changed if input structure changes
            hom_rev = nt_seq[int(modification[4]):int(modification[4])+30]
            hom_rev = hom_rev.reverse_complement()
            # Append the primer-bind overlaps
            proswp_seq = self.promoters[modification[2]]
            proswp_seq_revnoncom = str()  # storage string for the non-revcom'd of the back side of the PROSWP
            # Logic for processing the proswp sequence of capital letters to find the binding sites
            switch_trigger = True
            for nt in proswp_seq:
                if nt.isupper():
                    if switch_trigger:
                        hom_fwd += nt
                    else:
                        proswp_seq_revnoncom += nt
                else:
                    switch_trigger = False
            # revcom-ing the back end sequence of the proswp and adding them to the output list
            hom_rev = Seq(proswp_seq_revnoncom).reverse_complement() + hom_rev
            output_list.append(("hom_rev", hom_rev))
            output_list.append(("hom_fwd", hom_fwd))
        return output_list

A = AutoRefactor("/Users/bmendoza/Documents/Auto_Ref_Test/","PROSWPs")