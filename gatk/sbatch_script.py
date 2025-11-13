#!/usr/bin/env python
# coding: utf-8
import sys
import os
import argparse
import subprocess
import time

ENV = '/work/home/acfaa2ssz7/.bashrc'

def generate_sbatch_file(patient_id, job_name, job_out, job_err, job_mem, job_N, job_cpu, job_script,env='gatk4')->str:
    """
    Generate a sbatch file for submitting a job to the cluster.
    """
    
    sbatch_file_mem = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -o {job_out}
#SBATCH -e {job_err}
#SBATCH -p {patient_id}
#SBATCH -N {job_N}
#SBATCH -n {job_cpu}
#SBATCH --mem={job_mem}

module purge

source {ENV}
conda activate {env}

{job_script}
wait
"""
    sbatch_file = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -o {job_out}
#SBATCH -e {job_err}
#SBATCH -p {patient_id}
#SBATCH -N {job_N}
#SBATCH -n {job_cpu}

module purge

source {ENV}

conda activate {env}

{job_script}
wait
"""
    if job_mem == "0":
        with open(f"{job_name}.sbatch", "w") as f:
            f.write(sbatch_file)
    else:
        with open(f"{job_name}.sbatch", "w") as f:
            f.write(sbatch_file_mem)
    return f"{job_name}.sbatch"

def pysh(cmd) -> bool:
    '''
    执行输入的shell命令,并返回执行状态bool
    '''
    # p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    # stdout, stderr = p.communicate()
    # return_code = p.wait()
    subprocess.run(cmd, shell=True)
    # if return_code != 0:
    #     print(f'Error: {cmd} {stderr}')
    #     sys.exit(1)
    return True

def main():
    parser = argparse.ArgumentParser(description='Sbatch script generator')
    parser.add_argument('-i', '--input', type=str, required=True, help='input script file(eg:T1001\tpython /path/to/script.py -i ...)')
    parser.add_argument('-p', '--patient_id', type=str, required=True, help='patient id')
    parser.add_argument('-m', '--job_mem', type=str, required=True, help='memory size')
    parser.add_argument('-N', '--job_N', type=str, required=True, help='number of nodes')
    parser.add_argument('-n', '--job_cpu', type=str, required=True, help='number of cpu')
    parser.add_argument('-e', '--env', type=str, default='gatk4', help='environment name (default: gatk4)')
    args = parser.parse_args()
    
    input_file = args.input
    patient_id = args.patient_id
    job_mem = args.job_mem
    job_N = args.job_N
    job_cpu = args.job_cpu
    env = args.env

    job_script = {}
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                tag = line.split("\t")[0]
                script = line.split("\t")[1]
                if tag not in job_script:
                    job_script[tag] = script
                else:
                    print("Error: Duplicate tag found in input file:" + tag)
                    exit(1)
    
    sbatch_file_list = []
    for tag in job_script.keys():
        script_shell = job_script[tag]
        job_out = f"{tag}.out.%j"
        job_err = f"{tag}.err.%j"
        print(f'Submitted job for tag: {tag}')
        sbatch_file = generate_sbatch_file(patient_id,tag, job_out, job_err, job_mem, job_N, job_cpu, script_shell, env)
        sbatch_file_list.append(sbatch_file)

    for sbatch_file in sbatch_file_list:
        cmd = f'sbatch {sbatch_file}'
        time.sleep(2)
        pysh(cmd)
        
    
    print("All jobs submitted successfully.")


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'Error occurred: {e}')
        exit(1)

