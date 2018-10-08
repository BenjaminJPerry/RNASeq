import os
import collections


def read_in_tsv_file(infile):  # Make a function to parse the input files ###
    # type: (IO Stream) -> list


    """" (infile -> list[][])
        Reads in a tab seperated file, removes whitespace and " characters.
        Lines in infile -> rows in list[][]
        Tab delimited elements in infile -> columns in list[][]
    """
    input = infile
    data = []
    for line in input:
        items = line.strip().split('\t')  # Split the line into a string using tabs. Stripping.
        items = list([s.strip() for s in items])  # Stripping whitespace from each element
        items = list([s.strip('"') for s in items])
        data.append(items)
    return data


# Read in some ID table of VF39 genes with ID in column 1
RNAOld = 'VF39RNASeqIn.txt'
RNAIn = open(str(RNAOld), "Ur")
RNA_old = read_in_tsv_file(RNAIn)

# Read in VF39 ID:COG data
VF39COGIn = 'VF39COG.txt'
COGIn = open(str(VF39COGIn), "Ur")
COG = read_in_tsv_file(COGIn)

# ID Keyed dictionary of VF39 COGs
COG_dict = {}
COG_list = []
for entry in COG:  # make a locus_tag:entry2 dictionary for gff file
    if "#" in entry[0]:
        continue
    else:
        COG_dict.update({entry[0]: entry[9:11]})
        COG_list.append(entry[9])

# Update the RNASeq data file
RNANew = []
RNACOG_list = []

for gene in RNA_old:
    if '#' in gene[0]:
        continue
    elif gene[0] in COG_dict:
        new_entry = gene
        for i in COG_dict[gene[0]]:
            new_entry.append(i)
        RNANew.append(new_entry)
    else:
        continue

for gene1 in RNANew:
    RNACOG_list.append(gene1[7])

COG_hist = collections.Counter(COG_list)
RNACOG_hist = collections.Counter(RNACOG_list)
# make an output file
wkdir = os.getcwd()
outfile = str(wkdir + "VF39RNASeqData.updated.txt")
COG_out = open(outfile, "w")
print('Printing to ' + str(COG_out) + '\n')

# gff_new contains the updated gff file in a list[][]
for thing in RNANew:
    print_line = str("\t".join(thing))
    COG_out.write(print_line + '\n')


COG_out.write("### COG Genome Histogram ###\n")
for tag, count in COG_hist.items():
    COG_out.write('{}\t{}\n'.format(tag, count))

COG_out.write("### COG RNASeq Histogram ###\n")
for tag, count in RNACOG_hist.items():
    COG_out.write('{}\t{}\n'.format(tag, count))
COG_out.close()

exit()
