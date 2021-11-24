"""Code for identifying ideal capture sequences to generate hom arms and sgRNA sites.
SCRAPE: Simultaneous Cluster Recovery and Plasmid Exchange
INPUTS:
-Cluster to be pulled down as a contig ZID
-Vector to be assembled into
-BGC barriers (optional)"""

class SCRAPE:

    def __init__(self):
        self.cosmidDNA = str()
        self.borders = (0, 0)
        
