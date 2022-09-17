#! /usr/bin/env python3
''''
i will use this script to rename the transcript id in the gff file
chrUn   EVM     gene    2750227 2751450 .       +       .       ID=AS28G000180
chrUn   EVM     mRNA    2750227 2751450 .       +       .       ID=AS28G000180.1;Parent=AS28G000180
chrUn   EVM     exon    2750227 2751450 .       +       .       ID=AS28G000180.exon1;Parent=AS28G000180.1
chrUn   EVM     CDS     2750227 2751450 .       +       0       ID=cds.AS28G000180;Parent=AS28G000180.1
chrUn   EVM     CDS     2792358 2793425 .       +       0       ID=cds.AS28G000190;Parent=AS28G000190.1
chrUn   EVM     exon    2792358 2793425 .       +       .       ID=AS28G000190.exon1;Parent=AS28G000190.1
chrUn   EVM     gene    2792358 2794385 .       +       .       ID=AS28G000190
chrUn   EVM     mRNA    2792358 2794385 .       +       .       ID=AS28G000190.1;Parent=AS28G000190
and:
chr1A   EVM     gene    369603  371087  .       -       .       ID=AS01G000010;Name=AS01G000010
chr1A   EVM     mRNA    369603  371087  .       -       .       ID=AS01G000010;Parent=AS01G000010;Name=AS01G000010
chr1A   EVM     exon    369603  371087  .       -       .       ID=AS01G000010.exon1;Parent=AS01G000010
chr1A   EVM     CDS     369603  371087  .       -       0       ID=cds.AS01G000010;Parent=AS01G000010

order of the columns
0: chr
1: start
2: type
and rename the transcript id and the cds exon id
cds name and transcript name default
'''
import sys
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from tqdm import tqdm

# define the order of the categories
cat_size_order = CategoricalDtype(
    ['gene', 'mRNA', 'exon', 'CDS'],
    ordered=True
)


def main():
    input_gff: str = "A_strigosa_nc_sorted.gff3"
    output_gff_list: list = []
    file_line: int = 0
    input_gff_df = pd.read_csv(
        input_gff, sep="\t", header=None)
    # convert the 3rd column to categorical
    input_gff_df[2] = input_gff_df[2].astype(cat_size_order)
    input_gff_df.sort_values(by=[0, 3, 2], ascending=[
                             True, True, True], inplace=True)
    input_gff_list = np.array(input_gff_df)  # np.ndarray()
    input_gff_list = input_gff_list.tolist()  # list
    for line in tqdm(input_gff_list, desc="Processing annotation file"):
        if line[2] == "gene":
            # print(line)
            geneid = line[8].split(";")[0].split("=")[1]
            geneid_line = file_line
            output_gff_list.append(line)
        elif line[2] == "mRNA":
            parentid = line[8].split("Parent=")[1].split(";")[0]
            if geneid == parentid:
                # print("geneid: ", file_line, "parentid: ", geneid_line)
                transcript_number = file_line - geneid_line
                # print(transcript_number)
                if transcript_number == 1:
                    transcriptid = geneid + ".1"
                    # print(transcriptid)
                else:
                    print("Error: More than one transcript per gene")
                line[8] = "ID=" + transcriptid + ";Parent=" + \
                    geneid + ";Name=" + transcriptid
                output_gff_list.append(line)
            else:
                print("Error: mRNA does not match gene")
        elif line[2] == "exon":
            parentid = line[8].split("Parent=")[1]
            if geneid == parentid:
                line[8] = line[8].replace(
                    "Parent=" + parentid, "Parent=" + transcriptid)
                exon_name = line[8].split("ID=")[1].split(";")[
                    0].split(".")[1]
                line[8] = line[8].replace(
                    "ID=" + geneid + '.' + exon_name, "ID=" + transcriptid + "." + exon_name)
                output_gff_list.append(line)
            elif parentid == transcriptid:
                exon_name = line[8].split("ID=")[1].split(";")[
                    0].split(".")[1]
                line[8] = line[8].replace(
                    "ID=" + geneid + '.' + exon_name, "ID=" + transcriptid + "." + exon_name)
                output_gff_list.append(line)
            else:
                print('geneid: ', geneid, 'parentid: ',
                      parentid, 'transcriptid: ', transcriptid)
                print("Error: exon does not match gene")
        elif line[2] == "CDS":
            parentid = line[8].split("Parent=")[1]
            if geneid == parentid:
                line[8] = line[8].replace(
                    "Parent=" + parentid, "Parent=" + transcriptid)
                cds_name = line[8].split("ID=")[1].split(";")[
                    0].split(".")[0]
                line[8] = line[8].replace(
                    "ID=" + cds_name + '.' + geneid, "ID=" + transcriptid + "." + exon_name.replace("exon", "cds"))
                output_gff_list.append(line)
            elif parentid == transcriptid:
                cds_name = line[8].split("ID=")[1].split(";")[
                    0].split(".")[0]
                line[8] = line[8].replace(
                    "ID=" + cds_name + '.' + geneid, "ID=" + cds_name + "." + transcriptid)
                output_gff_list.append(line)
            else:
                print('geneid: ', geneid, 'parentid: ',
                      parentid, 'transcriptid: ', transcriptid)
                print("Error: cds does not match gene")
        else:
            print("Error: Unknown feature type")
        file_line += 1
    with open("A_strigosa_nc_sorted_transcript.gff3", "w") as output_gff:
        for line in output_gff_list:
            line[3] = str(line[3])
            line[4] = str(line[4])
            output_gff.write("\t".join(line) + "\n")
    # with open(input_gff, "r") as gff:
    #     for line in tqdm(gff, desc="Reading GFF"):
    #         if line.startswith("#"):
    #             continue
    #         else:


if __name__ == '__main__':
    #    import sys
    #    import os
    #    import re
    #    import argparse
    #
    #    parser = argparse.ArgumentParser(description='Rename transcript files')
    #    parser.add_argument('path', help='Path to transcript files')
    #    parser.add_argument('prefix', help='Prefix for transcript files')
    #    args = parser.parse_args()
    #
    #    path = args.path
    #    prefix = args.prefix
    #
    #    if not os.path.isdir(path):
    #        print('Path is not a directory')
    #        sys.exit(1)
    #
    #    files = os.listdir(path)
    #    files = [f for f in files if re.match(r'.*\.txt', f)]
    #
    #    for f in files:
    #        new_name = prefix + '_' + f
    #        os.rename(os.path.join(path, f), os.path.join(path, new_name))
    #
    #    print('Done')
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted')
        sys.exit(1)
