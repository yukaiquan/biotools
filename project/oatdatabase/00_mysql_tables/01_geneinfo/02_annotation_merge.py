#!/use/bin/env python
'''
author: yukaiquan
email: 1962568272@qq.com
usage: python 02_annotation_merge.py gff.file genome_mysql.file kog.file kegg.file nr.file uniprot.file output.file 
'''
import numpy as np
import sys
import datetime


def main(inputgffile: str, genomeid: str, kog_input: str, kegg_input: str, nr_input: str, uniport_input: str, outputfile: str):
    # debug
    # inputgffile = 'gff/A_strigosa_nc_sorted.gff3'
    # genomeid = '5'
    # kog_input = 'KOG/20220917output/strigosa.fasta_kog_outfmt6.txt_result'
    # kegg_input = 'KEGG/20220915/strigosa.fasta.querry2KO.gz.koid'
    # nr_input = 'nr/01_annotation_0916/blast_annotation/strigosa.fasta.out.ann'
    # uniport_input = 'uniprot/20220919output/strigosa.fastaoutfmt6.txt_anndetail.tsv'
    # outputfile = 'test_out.txt'
    '''
    inputgffile: input gff file
    format: seqid, source, type, start, end, score, strand, phase, attributes
    chr1A	EVM	gene	251998	253535	.	-	.	ID=A.satnudSFS1A01G007012
    chr1A	EVM	mRNA	251998	253535	.	-	.	ID=A.satnudSFS1A01G007012.1;Parent=A.satnudSFS1A01G007012
    chr1A	EVM	CDS	251998	252744	.	-	0	ID=CDS2.A.satnudSFS1A01G007012.1;Parent=A.satnudSFS1A01G007012.1
    chr1A	EVM	exon	251998	252744	.	-	.	ID=exon2.A.satnudSFS1A01G007012.1;Parent=A.satnudSFS1A01G007012.1
    chr1A	EVM	CDS	252870	253535	.	-	0	ID=CDS1.A.satnudSFS1A01G007012.1;Parent=A.satnudSFS1A01G007012.1
    chr1A	EVM	exon	252870	253535	.	-	.	ID=exon1.A.satnudSFS1A01G007012.1;Parent=A.satnudSFS1A01G007012.1
    kog_input: input kog file
    format: seqid, kog
    A.satnudSFSUN01G001420	At1g35420	[R] KOG3043 Predicted hydrolase related to dienelactone hydrolase
    A.satnudSFSUN01G001416	At4g23240	[T] KOG1187 Serine/threonine protein kinase
    A.satnudSFSUN01G001415	At4g23240	[T] KOG1187 Serine/threonine protein kinase
    A.satnudSFSUN01G001414	At4g23240	[T] KOG1187 Serine/threonine protein kinase
    kegg_input: input kegg file
    format: seqid, kegg
    A.satnudSFSUN01G001420 K01061   99.67  118.4   2.9e-33 carboxymethylenebutenolidase [EC:3.1.1.45]
    A.satnudSFSUN01G001419 K09833  289.43  450.9  3.3e-134 homogentisate phytyltransferase / homogentisate geranylgeranyltransferase [EC:2.5.1.115 2.5.1.116]
    A.satnudSFSUN01G001413 K08176  337.33  856.3  1.1e-256 MFS transporter, PHS family, inorganic phosphate transporter
    nr_input: input nr file
    format: seqid, nr
    A.satnudSFS1A01G000002	A0A8I7B245	69.8	440	7.3e-169	600.1
    A.satnudSFS1A01G000005	A0A8I6WNZ6	63.9	765	3.2e-231	807.7
    A.satnudSFS1A01G000009	A0A3B5Y8V4	93.8	433	3.7e-232	810.1
    uniport_input: input uniport file
    format: seqid, uniport
    A.satnudSFS1A01G000002	KAI4968384.1	70.0	440	4.5e-169	600.9	hypothetical protein ZWY2020_058039 [Hordeum vulgare]
    A.satnudSFS1A01G000005	QDB64260.1	95.0	756	0.0e+00	1153.7	early flowering 3 [Avena sativa]
    A.satnudSFS1A01G000007	XP_047065255.1	60.8	222	8.5e-65	254.2	increased DNA methylation 1-like [Lolium rigidum]
    '''
    # Open the input file
    # gff_ndarray = np.genfromtxt(
    #     inputgfffile, delimiter="\t", dtype=None, names=True)
    # gff_ndarray[gff_ndarray[2] == "gene"]
    gff_gene_list: list = []
    # Open the input file
    with open(inputgffile, "r") as gff_file:
        # Read the file line by line
        for line in gff_file:
            if line.startswith("#"):
                continue
            # Split the line into a list
            line_list = line.split("\t")
            # Check if the line is a gene
            if line_list[2] == "gene":
                if ';' in line_list[8]:
                    line_list[8] = line_list[8].strip().split(';')[
                        0].split('=')[1]
                    gff_gene_list.append(line_list)
                else:
                    line_list[8] = line_list[8].strip().split('=')[1]
                    gff_gene_list.append(line_list)

    gff_gene_list = np.array(gff_gene_list)
    # replace the gene name with the gene id
    # gff_gene_list[:, 8] = np.array([geneid[0].split(
    #     '=')[1] for geneid in np.char.split(gff_gene_list[:, 8], sep=";")])
    # gff_gene_list[:, 8] = np.char.replace(gff_gene_list[:, 8], "ID=", "")
    # gff_gene_list[:, 8] = np.char.replace(gff_gene_list[:, 8], "\n", "")
    # caculate gene length
    gff_gene_list = np.append(gff_gene_list, abs(gff_gene_list[:, 4].astype(
        float) - gff_gene_list[:, 3].astype(float)).reshape(-1, 1), axis=1)
    # replace +- with 1 and - with 0
    gff_gene_list[:, 6] = np.where(gff_gene_list[:, 6] == "+", 1, 0)
    # add a column for create time and update time
    arry_len = len(gff_gene_list)
    gff_gene_list = np.append(gff_gene_list, np.char.replace(
        np.zeros((arry_len, 1), dtype=str), '', '2022-09-27 23:03:00.000000'), axis=1)
    gff_gene_list = np.append(gff_gene_list, np.char.replace(
        np.zeros((arry_len, 1), dtype=str), '', '2022-09-27 23:03:00.000000'), axis=1)

    gff_gene_list = np.delete(gff_gene_list, [1, 2, 5, 7], axis=1)

    kog_dict = {}
    # read KOG result
    with open(kog_input, "r") as kog_file:
        kog_list = kog_file.readlines()
        for line in kog_list:
            line_list = line.split("\t")
            kog_ann = line_list[2].split()[1]
            kog_dict[line_list[0]] = kog_ann
    del kog_list

    # add KOG annotation to gff_gene_list
    kog_arry = np.array([kog_dict.get(gene_id)
                        for gene_id in gff_gene_list[:, 4]])
    gff_gene_list = np.append(gff_gene_list, kog_arry.reshape(-1, 1), axis=1)
    del kog_arry

    # KEGG annotation
    kegg_dict = {}
    # read KEGG result
    with open(kegg_input, "r") as kegg_file:
        kegg_list = kegg_file.readlines()
        for line in kegg_list:
            line_list = line.strip().split("\t")
            # print(line_list)
            kegg_dict[line_list[0]] = line_list[1]

    kegg_arry = np.array([kegg_dict.get(gene_id)
                          for gene_id in gff_gene_list[:, 4]])
    gff_gene_list = np.append(gff_gene_list, kegg_arry.reshape(-1, 1), axis=1)
    del kegg_arry

    # NR annotation
    nr_dict = {}
    # read NR result
    with open(nr_input, 'r') as f:
        nr_list = f.readlines()
        for line in nr_list:
            line_list = line.split("\t")
            nr_dict[line_list[0]] = line_list[1]

    nr_arry = np.array([nr_dict.get(gene_id)
                       for gene_id in gff_gene_list[:, 4]])
    gff_gene_list = np.append(gff_gene_list, nr_arry.reshape(-1, 1), axis=1)
    del nr_arry

    # uniport annotation
    uniport_dict = {}
    with open(uniport_input, 'r') as f:
        uniport_list = f.readlines()
        for line in uniport_list:
            line_list = line.split("\t")
            uniport_dict[line_list[0]] = line_list[1]

    uniport_arry = np.array([uniport_dict.get(gene_id)
                            for gene_id in gff_gene_list[:, 4]])
    gff_gene_list = np.append(
        gff_gene_list, uniport_arry.reshape(-1, 1), axis=1)

    # genome id add
    genome_id = np.char.replace(
        np.zeros((arry_len, 1), dtype=str), '', genomeid)
    gff_gene_list = np.append(gff_gene_list, genome_id, axis=1)

    # gene description
    gene_desc = np.char.replace(
        np.zeros((arry_len, 1), dtype=str), '', 'gene first description')
    gff_gene_list = np.append(gff_gene_list, gene_desc, axis=1)
    gff_gene_list = gff_gene_list[:, [
        4, 0, 1, 2, 3, 5, 13, 6, 7, 12, 9, 8, 10, 11]]
    gff_gene_list = gff_gene_list[np.argsort(gff_gene_list[:, 0])]

    # write to csv
    np.savetxt(outputfile, gff_gene_list, delimiter=",", fmt="%s")


if __name__ == "__main__":
    try:
        input_gff = sys.argv[1]
        genome_mysql_id = str(sys.argv[2])
        kog_result = sys.argv[3]
        kegg_result = sys.argv[4]
        nr_result = sys.argv[5]
        uniport_result = sys.argv[6]
        output_csv = str(sys.argv[7]).strip()
        print('#'*50)
        start_time = datetime.datetime.now()
        print('start_time: ', start_time)
        print('input_gff: ', input_gff, 'genome_mysql_id: ', genome_mysql_id, 'kog_result: ', kog_result, 'kegg_result: ',
              kegg_result, 'nr_result: ', nr_result, 'uniport_result: ', uniport_result, 'output_csv: ', output_csv)
        main(inputgffile=input_gff, genomeid=genome_mysql_id, kog_input=kog_result,
             kegg_input=kegg_result, nr_input=nr_result, uniport_input=uniport_result, outputfile=output_csv)
        print('end_time: ', datetime.datetime.now())
        print('cost_time: ', datetime.datetime.now() - start_time)
        print('#'*50)
    except Exception as e:
        print(e)
        sys.exit(1)
