import sys
import gzip

with gzip.open(sys.argv[1], "rb") as f:
    for line in f:
        lsplit = line.decode('utf-8').rstrip().split("\t")
        if lsplit[7] and lsplit[0]:
            line = lsplit[0] + "\t" + lsplit[7]
            print(line)
