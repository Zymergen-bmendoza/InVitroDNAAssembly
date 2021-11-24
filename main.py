# auto design the gibson assembler to generate the overhangs with the optimal GC content.
# auto design sgRNAs to cut at the appropriate spot, minimizing issues of off-target
# for multiplexing: how big of pieces are we dealing with?

"""Quick psuedocode for reference:
Get genbank file, download"""

import os
from Bio import SeqIO
import re

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

        # Import the clusters via annotation parameters
        #self.parse_clusters(workingdir)

        # Now go through the clusters that have been imported:
        output_dict = dict()
        for cluster in self.params_dict:
            output_dict[cluster] = list()
            output_dict[cluster].append(self.guide_selector(cluster))
            if self.mod_style == "PROSWP":
                output_dict[cluster].append(self.hifi_hom_arm_generator(cluster))
                self.output(output_dict)
            elif self.mod_style == "BBSWP":
                output_dict[cluster].append(self.hifi_hom_arm_generator(cluster))
                self.output(output_dict)
            elif self.mod_style == "GENEDEL":
                output_dict[cluster].append(self.hifi_hom_arm_generator(cluster))
                self.output(output_dict)
            else:
                print("Invalid Modification selected.")


    """Function opens the output stream to the file for building the factory order."""
    def output(self, out):
        out_file_name = "myoutput_" + self.mod_style + ".csv"
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

        for line in f:
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
        y = open("promoter_parts.csv")
        for line in y:
            prom = line[:-1].split("\t")
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

    """Function for once the location is established, the guides can be designed.
        ideal_location input is where the cut site (3 bp away from NGG) should be."""
    def guide_selector(self, selected_cluster,):
        # Getting the information from the parameter dictionary:
        cluster_params = self.params_dict[selected_cluster]

        # Output list:
        guide_list = list()

        locations = list()
        for refsite in cluster_params:
            locations.append(int(refsite[2]))
            locations.append(int(refsite[3]))

        for ideal_location in locations:
            nt_seq = SeqIO.read(selected_cluster + ".gb", "genbank").seq
            if self.mod_style == "BBSWP":
                gg_seq = str(nt_seq[ideal_location - 30:ideal_location + 30])
                cc_seq = str(nt_seq[ideal_location - 30:ideal_location + 30])
            else:
                # for PROSWP/GENEDEL ideal location math for GG: + 14 or - 6 for CC it is the opposite (- 14 or + 6)
                gg_seq = str(nt_seq[ideal_location - 6:ideal_location + 14])
                cc_seq = str(nt_seq[ideal_location - 14:ideal_location + 6])

            # Now add them to lists to parse through
            locs_c = [m.start() for m in re.finditer('CC', cc_seq)]
            locs_g = [m.start() for m in re.finditer('GG', gg_seq)]

            # Now correct for the cut site and find the best one in each of the 2 lists
            ideal_guide_diff = 100  # arbitrarily large number to automatically assign the first guide
            ideal_guide_loc = int()
            revcom_flag = False
            # For NGG Pams:
            for loc in locs_g:
                diff = abs(loc-10)  # corrected for start of subsequence and cut site
                if diff < ideal_guide_diff:
                    ideal_guide_loc = loc
                    ideal_guide_diff = diff
            # For revcom Pams:
            for loc in locs_c:
                diff = abs(loc-10)  # corrected for start of subsequence and cut site
                if diff < ideal_guide_diff:
                    ideal_guide_loc = loc
                    ideal_guide_diff = diff
                    revcom_flag = True

            # Now get the ideal guide and get its true sequence from nt_seq
            # Math needed for this: need to get the ideal location and set it against the actual location then backtrack
            if not revcom_flag:
                srtloc = ideal_location - 6
                guide_list.append((nt_seq[srtloc+ideal_guide_loc-21:srtloc+ideal_guide_loc-1],
                                   srtloc+ideal_guide_loc-4))
            else:
                srtloc = ideal_location - 14
                guide_list.append((nt_seq[srtloc + ideal_guide_loc+3:srtloc + ideal_guide_loc+23],
                                   srtloc+ideal_guide_loc+6))

        return guide_list

    """Inputs: The parameters and parts files, the genbank of the cluster.
        Output: The primer sequences for generating the homology arms."""
    def hifi_hom_arm_generator(self, cluster):
        output_list = list()
        nt_seq = SeqIO.read(cluster + ".gb", "genbank").seq
        for proswp in self.params_dict[cluster]:
            hom_rev = nt_seq[int(proswp[2])-30:int(proswp[2])]
            hom_rev = hom_rev.reverse_complement()
            hom_fwd = nt_seq[int(proswp[3]):int(proswp[3])+30]
            output_list.append(("hom_rev", hom_rev))
            output_list.append(("hom_fwd", hom_fwd))
        return output_list

    def bbswap_primers_generator(self, cluster):
        output_list = list()
        nt_seq = SeqIO.read(cluster + ".gb", "genbank").seq
        for bbswp_cut in self.params_dict[cluster]:
            hom_pcr_rev = nt_seq[int(bbswp_cut[3]):int(bbswp_cut[3])+30]
            hom_pcr_rev = hom_pcr_rev.reverse_complement()
            hom_pcr_fwd = nt_seq[int(bbswp_cut[2])-30:int(bbswp_cut[2])]
            output_list.append(("hompcrrev", hom_pcr_rev))
            output_list.append(("hompcrfwd", hom_pcr_fwd))
        return output_list



# Options for style: PROSWP, BBSWP, GENEDEL
A = AutoRefactor("C:/Users/bmendoza/Documents/Auto_Ref_Test/", mod_style="PROSWP")


        


