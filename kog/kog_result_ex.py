#! /usr/bin/env python3
'''
use:
python kog_result_ex.py kog kog_result.txt kog_result_ex.txt
kog_result.txt:
A.satnudSFSUN01G001420  At1g35420       48.1    243     83      6       129     333     73      310     1.1e-53 209.1
A.satnudSFSUN01G001416  At4g23240       45.5    343     154     5       239     554     11      347     5.3e-75 280.8
A.satnudSFSUN01G001416  At4g00960       43.5    338     174     8       243     572     43      371     1.9e-64 245.7
A.satnudSFSUN01G001416  At1g70740       41.0    300     169     5       244     539     6       301     4.2e-56 218.0
A.satnudSFSUN01G001416  At2g20300       37.2    285     171     4       244     525     28      307     5.9e-50 197.6
kog_result_ex.txt:
A.satnudSFSUN01G001420  At1g35420       [R] KOG3043 Predicted hydrolase related to dienelactone hydrolase
A.satnudSFSUN01G001416  At4g23240       [T] KOG1187 Serine/threonine protein kinase
A.satnudSFSUN01G001415  At4g23240       [T] KOG1187 Serine/threonine protein kinase
A.satnudSFSUN01G001414  At4g23240       [T] KOG1187 Serine/threonine protein kinase
A.satnudSFSUN01G001413  At2g38940       [P] KOG0252 Inorganic phosphate transporter
A.satnudSFSUN01G001412  At1g17250       [R] KOG0619 FOG: Leucine rich repeat
A.satnudSFSUN01G001411  At1g78530       [T] KOG1187 Serine/threonine protein kinase
A.satnudSFSUN01G001409  At1g67580       [R] KOG0663 Protein kinase PITSLRE and related kinases
'''
import re
import sys
from tqdm import tqdm

# input_kog_re = 'ACD.SFS.pep.fasta_kog_outfmt6.txt'
# output_kog_re = 'ACD.SFS.pep.fasta_kog_outfmt6_result.txt'
kog_class = sys.argv[1]
input_kog_re = sys.argv[2]
output_kog_re = sys.argv[3]

id_kog = {}
kog_list = []
p = re.compile(r'^$')

with open(kog_class, 'r') as f:
    for line in tqdm(f, desc='kog_class read'):
        if re.findall(p, line):
            continue
        else:
            kog_list.append(line.strip())

for line in tqdm(kog_list, desc='kog_class to dict'):
    if line.startswith('['):
        kog_ann = line.strip()
    else:
        # print(line)
        try:
            id_kog[line.split(':')[1].strip()] = kog_ann
        except IndexError:
            print('index error content(but dont care) %s' % line)
del kog_list

kog_result = {}
with open(input_kog_re, 'rb') as f:
    for line in tqdm(f, desc='kog_result read'):
        line = line.decode('utf-8')
        if line.startswith('#'):
            continue
        else:
            line = line.strip().split('\t')
            if float(line[10]) <= 1e-5:
                if line[0] not in kog_result:
                    kog_result[line[0]] = [line[1], id_kog[line[1]]]
                else:
                    continue
with open(output_kog_re, 'w') as f:
    for key, value in tqdm(kog_result.items(), desc='kog_result write to file'):
        f.write('%s\t%s\t%s\n' % (key, value[0], value[1]))
