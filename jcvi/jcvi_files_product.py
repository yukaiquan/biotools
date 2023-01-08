#!/use/bin/env python3

import sys
from subprocess import Popen
import time


def path_remake(path):
    return path.replace(' ', '\ ').replace('(', '\(').replace(')', '\)')


def bed_product(gff_file, file_name) -> bool:
    cmd = f"python -m jcvi.formats.gff bed --type=mRNA --key=ID {gff_file} -o {file_name}.bed"
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return True


def bed_uniq(file_name) -> bool:
    cmd = f"python -m jcvi.formats.bed uniq {file_name}.bed"
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return True


def mv_bed(file_name) -> bool:
    cmd = f"mv {file_name}.uniq.bed {file_name}.bed"
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return True


def get_seq_id(file_name) -> bool:
    cmd = f"cut -f 4 {file_name}.bed > {file_name}.id"
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return True


def get_seq(seq_file, file_name, out_file) -> bool:
    cmd = f"seqkit grep -f {file_name}.id {path_remake(seq_file)} | seqkit seq -i > {out_file}"
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return True


def main(gff_file, cds_file, pep_file, file_name):
    start = time.time()
    if bed_product(gff_file, file_name):
        print("bed product done")
    else:
        print("bed product failed")
        sys.exit(1)
    if bed_uniq(file_name):
        print("bed uniq done")
    else:
        print("bed uniq failed")
        sys.exit(1)
    if mv_bed(file_name):
        print("mv bed done")
    else:
        print("mv bed failed")
        sys.exit(1)
    if get_seq_id(file_name):
        print("get id done")
    else:
        print("get id failed")
        sys.exit(1)
    if get_seq(cds_file, file_name, f"{file_name}.cds"):
        print("get cds done")
    else:
        print("get cds failed")
        sys.exit(1)
    if get_seq(pep_file, file_name, f"{file_name}.pep"):
        print("get pep done")
    else:
        print("get pep failed")
        sys.exit(1)
    end = time.time()
    print(f"Total time: {end - start} s")


if __name__ == "__main__":
    try:
        gff_file = sys.argv[1]
        cds_file = sys.argv[2]
        pep_file = sys.argv[3]
        file_name = sys.argv[4]
        main(gff_file, cds_file, pep_file, file_name)
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me!")
