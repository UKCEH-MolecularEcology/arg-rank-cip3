# MODULES
import os, re
import glob

cip3_samples = [line.strip() for line in open("cip3_samples").readlines()]

rule all:
    input:
        expand(os.path.join("output/{sample}/{sample}_out.txt"), sample = cip3_samples)

rule mc:
    input:
        "rawdata/{sample}.fna"
    output:
        "output/{sample}/{sample}_out.txt"
    log:
        "logs/{sample}_mc.log"
    conda:
        "/home/roseelliott/miniforge3/envs/metacompare"
    shell:
        "mkdir -p $(dirname {output}) && "
        "python metacompare.py -c {input} -o $(dirname {output})"

