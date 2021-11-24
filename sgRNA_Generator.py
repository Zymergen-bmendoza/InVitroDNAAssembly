"""File for generating Hitpick plans for sgRNA synthesis."""

import os
from Bio import SeqIO
import re

"""Input: Modification style is set by specific class that is calling the function.  Everytime guide_selector function
is called, a new cluster site must be selected
Output: A list of guides for the cluster site based on the type of modification."""


class sgRNAGenerator:

    def __init__(self, mod_style, filedir):
        self.mod_style = mod_style
        os.chdir(filedir + mod_style)

    """Function for once the location is established, the guides can be designed.
            ideal_location input is where the cut site (3 bp away from NGG) should be."""

    def guide_selector(self, selected_cluster, locations):
        # Output list:
        guide_list = list()

        for ideal_loc in locations[-2:]:  # -2 to get us the last two entries which will be the cut sites
            # change ideal location to numeric:
            ideal_location = int(ideal_loc)
            nt_seq = SeqIO.read(selected_cluster + ".gb", "genbank").seq
            if ideal_location < 0:
                ideal_location = len(nt_seq) + ideal_location
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
                diff = abs(loc - 10)  # corrected for start of subsequence and cut site
                if diff < ideal_guide_diff:
                    ideal_guide_loc = loc
                    ideal_guide_diff = diff
            # For revcom Pams:
            for loc in locs_c:
                diff = abs(loc - 10)  # corrected for start of subsequence and cut site
                if diff < ideal_guide_diff:
                    ideal_guide_loc = loc
                    ideal_guide_diff = diff
                    revcom_flag = True

            # Now get the ideal guide and get its true sequence from nt_seq
            # Math needed for this: need to get the ideal location and set it against the actual location then backtrack
            if self.mod_style == "BBSWP":
                if not revcom_flag:
                    srtloc = ideal_location - 30
                    guide_list.append((nt_seq[srtloc + ideal_guide_loc - 21:srtloc + ideal_guide_loc - 1],
                                       srtloc + ideal_guide_loc - 4))
                else:
                    srtloc = ideal_location - 30
                    guide_list.append((nt_seq[srtloc + ideal_guide_loc + 3:srtloc + ideal_guide_loc + 23].reverse_complement(),
                                       srtloc + ideal_guide_loc + 6))
            else:
                if not revcom_flag:
                    srtloc = ideal_location - 6
                    guide_list.append((nt_seq[srtloc + ideal_guide_loc - 21:srtloc + ideal_guide_loc - 1],
                                       srtloc + ideal_guide_loc - 4))
                else:
                    srtloc = ideal_location - 14
                    guide_list.append((nt_seq[srtloc + ideal_guide_loc + 3:srtloc + ideal_guide_loc + 23].reverse_complement(),
                                       srtloc + ideal_guide_loc + 6))

        return guide_list
