
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv, os


"""Function takes the sequences from a .csv file and turns them into genbank files in a directory."""
def sgRNAs_to_gb(mydir, sgrna_file):
    os.chdir(mydir)
    with open(sgrna_file, 'r') as f:
        guides = csv.DictReader(f, delimiter="\t")
        for item in guides:
            print(item)
            # Create a sequence
            sequence_string = item['Target Sequence'] +\
                              "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
            if not sequence_string.startswith('G'):
                sequence_string = "G" + sequence_string
            sequence_object = Seq(sequence_string)


            # Create a record
            record = SeqRecord(sequence_object,
                               id=item['sgRNAID'], # random accession number
                               name=item['Unique Name'],
                               description=item['Oligo Name'],
                               annotations={"molecule_type": "DNA"})

            # Add annotation
            feature = SeqFeature(FeatureLocation(start=1, end=20), type='misc_RNA')
            record.features.append(feature)

            # Save as GenBank file
            name = item['Unique Name'] + ".gb"
            output_file = open(name, 'w')
            SeqIO.write(record, output_file, 'genbank')
            output_file.close()


def smallDNAcomponent(path, parts_list):
    os.chdir(path)
    with open(parts_list, 'r') as f:
        guides = csv.DictReader(f, delimiter="\t")
        for item in guides:
            print(item)
            # Create a sequence
            sequence_object = Seq(item['Sequence'])

            # Create a record
            record = SeqRecord(sequence_object,
                               id=item['DNApartID'],  # random accession number
                               name=item['Part Name'],
                               description=item['Part Description'],
                               annotations={"molecule_type": "DNA"})

            # Save as GenBank file
            name = item['Unique Name'] + ".gb"
            output_file = open(name, 'w')
            SeqIO.write(record, output_file, 'genbank')
            output_file.close()


#sgRNAs_to_gb("/Users/bmendoza/Documents/sgRNAdatabase/","sgRNALIMSlist.tsv")
smallDNAcomponent("/Users/bmendoza/Documents/DNAcomponents/","uploadlist.tsv")