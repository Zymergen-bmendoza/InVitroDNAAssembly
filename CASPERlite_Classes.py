"""This file, imported into CASPERlite, stores the gRNA and PAMeval classes necessary for execution."""

import re
import os


class CASPERlite:

    def __init__(self, file_path, fa_file, scr_cutoff):
        self.fa_sequences = dict()  # for storing the locations of the cluster .fa files
        f = open(file_path + fa_file)
        # Code for adding the fasta sequences to the batch.
        s_id = str()
        myseq = ""
        for line in f:
            if line.startswith(">"):
                if len(myseq) > 1:
                    self.fa_sequences[s_id] = myseq
                    myseq = ""
                    s_id = line[:-1]
                else:  # for the first line exception to the pattern
                    s_id = line[:-1]
            else:
                myseq += line[:-1]
        # for adding last cluster:
        self.fa_sequences[s_id] = myseq
        f.close()

        #  CRISPRscan data import
        self.cscan = list()
        self.import_crisprscan()

        #  Stores the targets for output
        self.targets = list()  # stores the target sequences
        self.scr_cutoff = scr_cutoff

        # Overall run function
        self.run(file_path)

    def run(self, file_path):
        os.chdir(file_path)
        out = open(file_path + "OUTFILECH5acontigs" + ".txt", 'w')
        for cluster_seq in self.fa_sequences:
            print(cluster_seq)
            self.run_on_target_cl(self.fa_sequences[cluster_seq])
            # Sorts and prioritizes targets based on rules
            #sortedtargets = self.sort_and_prioritize(len(self.fa_sequences[cluster_seq], 30)
            sortedtargets = sorted(self.targets, key=lambda x: x[2], reverse=True)
            for i in range(25):
                out.write(cluster_seq + '\t')
                out.write(sortedtargets[i][0] + '\t')
                out.write(str(sortedtargets[i][1]) + '\t')
                out.write(str(sortedtargets[i][2]) + '\t')
                out.write(str(sortedtargets[i][3]) + '\n')
            self.targets.clear()
        out.close()

    #  This is the simplified version of CASPER on-target algorithm
    def run_on_target_cl(self, cluster):
        for match in re.finditer("GG", cluster[21:]):
            pam = match.start() + 21
            pass_target = self.trim_targets(cluster[pam-21:pam-1], pam, "+")
            # Error checking before adding to the target vector:
            if len(pass_target[0]) == 20 and pass_target[2] > self.scr_cutoff:
                self.targets.append(pass_target)
        # Now check the revcom
        rcluster = self.revcom(cluster)
        for match in re.finditer("GG", rcluster[21:]):
            pam = match.start() + 21
            pass_target = self.trim_targets(rcluster[pam - 21:pam - 1], pam, "-")
            # Error checking before adding to the target vector:
            if len(pass_target[0]) == 20 and pass_target[2] > self.scr_cutoff:
                self.targets.append(pass_target)

    #  This is the simplified version of the on-target sequence optimization from the CASPER algorithm
    def trim_targets(self, put_target, pam, direction):
        score = 50
        # Check to see if the target has a string of T's (this is bad for IVT)
        if put_target.find("TTTTT") != -1:
            return put_target, direction, 0, pam

        # Score with CRISPRscan right here
        for pattern in self.cscan:
            index = int(pattern[1]) - 1
            if index >= 20 or index < 0:
                continue
            if pattern[0][1] == 'x':
                if put_target[index] == pattern[0][0]:
                    score *= (1 + float(pattern[2])*2)
            else:
                if put_target[index:index+1] == pattern[0]:
                    score *= (1 + float(pattern[2])*2)

        # check to make sure there aren't too many PAMs already in the sequence
        factor_val = 0
        rev_factor_val = 0
        for match in re.finditer("GG", put_target):
            factor_val += 1
        score -= (10*factor_val)
        for match in re.finditer("CC", put_target):
            rev_factor_val += 1
        score -= (5*rev_factor_val)
        if put_target.startswith("G"):
            score *= 1.10

        # Confirm acceptable GC content (lower than 75%)
        if put_target.count("G") + put_target.count("C") > 14:
            score = score/2

        return put_target, direction, score, pam

    """def sort_and_prioritize(self, c_length, num_targets):
        ret_vec = list()
        sortedtargets = sorted(self.targets, key=lambda x: x[1], reverse=True)
        # Grab pairs for primers
        numspans = c_length/3000
        for sgrna in sortedtargets:
            if sgrna[2] """


    def import_crisprscan(self):
        f = open("C:/Users/bmendoza/Desktop/D1_CASPER/CASPERinfo.txt")
        inline = False
        for line in f:
            if inline:
                if line.startswith("-"):
                    break
                nucleotuple = tuple(line[:-1].split('\t'))
                self.cscan.append(nucleotuple)
            elif line.startswith("CRISPRSCAN"):
                inline = True

        f.close()

    def revcom(self, seq):
        rev = {"A": "T", "C": "G", "T": "A", "G": "C"}
        retseq = ""
        for base in seq:
            retseq = rev[base] + retseq
        return retseq


C = CASPERlite("C:/Users/bmendoza/Desktop/CASPERcontigs/", "CH52021contigs.fa", 50)
