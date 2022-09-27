import sys

# inputfile = 'ACD_sativa_pep.fasta'
inputfile = sys.argv[1]
fa_dict = {}
with open(inputfile, 'r') as f:
    for line in f:
        if line.startswith('>'):
            fa_id = line.strip()
            fa_dict[fa_id] = ''
        else:
            fa_dict[fa_id] += line.strip()

fa_name = []

for key in fa_dict:
    fa_name.append(key)

n = 45000
output = [fa_name[i:i + n] for i in range(0, len(fa_name), n)]


def write_file(output, fa_dict, fa_name):
    for i in range(len(output)):
        with open('%s_%s.fasta' % (inputfile, i), 'w') as f:
            for j in output[i]:
                f.write(j + '\n')
                f.write(fa_dict[j] + '\n')


write_file(output, fa_dict, fa_name)
