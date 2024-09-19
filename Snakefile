'''
Author: Txj, Yzy & Lsy
Affiliation: HumPOG, Fudan University
Date: Nov 23 2023
Aim:  Snakemake workflow for paired-end NGS reads mapping & variant-calling with CPC-graph.
#Requirements: vg==1.49, samtools, gatk, freebayes(conda), bedtools, abra2, deepvariant
Run:  snakemake  -s Snakefile
Latest modification: August 4th 2023
'''
import os

# Load configuration from config.yaml

configfile: "config.yaml"
filelist = config['filelist']
wd = config['work_dir']

# linear reference
ref = config['ref']
size_file = config['size_file']
ref_name = config['ref_name']

# PanGenome graph
graph_prefix=config['graph_prefix']
dist_file = config['dist_file']
min_file = config['min_file']
gbz_file = config['gbz_file']
# pangenie
Graph_preprocessed_vcf = config['Graph_preprocessed_vcf']

# enviroment
threads = config['threads']
run_deepvariant = config['run_deepvariant']
run_glnexus = config['run_glnexus']
run_pangenie = config['run_pangenie']

os.system('mkdir -p '+wd)

# Function to read in the filelist
from functions.basic import parse_filelist
# Function to get the chromsomes names from the graph reference
from functions.basic import get_chromosomes
# Function to split the reference
from functions.basic import split_chromosomes

# Process the filelist
try:
    fq_dict, pair_flag = parse_filelist(filelist)
except (FileNotFoundError, ValueError) as e:
    print(e)

samples = fq_dict.keys()


# Get the reference chromosomes
CHROMOSOMES = ['chr'+i for i in list(map(str,range(1,23)))+["X","Y","M"]]
os.popen('grep > '+ref+' | cut -d">" -f2')


with open(wd+'/paths.list','w') as f:
    f.writelines([graph_prefix+i+"\n" for i in CHROMOSOMES])




rule all:
    input:
        # final checkpoint for every sample
        expand(wd + "/logs/{sample}.checkpoints.end", sample = samples),

rule Split_ref:
    input:
        ref = ref,
        Graph_preprocessed_vcf = Graph_preprocessed_vcf
    output:
        ref_chr = wd + "/temp/ref/"+ref_name+".{chr}.fa",
        ref_chr_index = wd + "/temp/ref/"+ref_name+".{chr}.fa.fai",
        ref_vcf_chr = wd + "/temp/ref/"+ref_name+".{chr}.pangenie.vcf",
    log: wd + "/logs/00.split_ref.{chr}."+ref_name+".log"
    threads: 1
    resources: tmpdir= wd + "/temp/tmp", mem_mb=1024
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        samtools faidx {input.ref} {wildcards.chr} > {output.ref_chr}
        samtools faidx {output.ref_chr}
        bcftools view -r {wildcards.chr} {input.Graph_preprocessed_vcf} -o {output.ref_vcf_chr}
        """


#Step1: Run vg giraffe for graph reads mapping (~150 CPU hours, 5.5 hours with 36 threads) 
rule Graph_Mapping:
    input:
        fq = lambda wildcards: fq_dict[wildcards.sample],
    output:
        dist = wd + "/temp/{sample}/pangenome.dist",
        min = wd + "/temp/{sample}/pangenome.min",
        gbz = wd + "/temp/{sample}/pangenome.gbz",
        gam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam",
        gam_stat = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.gam_stat.txt"
    log: wd + "/logs/01.graph_mapping.{sample}."+ref_name+".log"
    benchmark: wd + "/benchmarks/01.graph_mapping.{sample}."+ref_name+".txt"
    threads: threads
    resources: tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*62
    run:
        shell("mkdir -p {wd}/temp/{wildcards.sample}/tmp")
        shell("cp {dist_file} {output.dist}")
        shell("cp {min_file} {output.min}")
        shell("cp {gbz_file} {output.gbz}")
        # vg giraffe
        cmd_template = "(vg giraffe -N {wildcards.sample} --read-group 'SM:{wildcards.sample}' -t {threads} -d {output.dist} -m {output.min} -Z {output.gbz} -f "
        fq_cmd = " -f ".join(input.fq)
        shell(cmd_template + fq_cmd + ' > {output.gam}) > {log} 2>&1')
        # vg stat
        cmd = "vg stats --threads 16 -a {output.gam} >> {output.gam_stat}"
        shell(cmd)

# ruleorder: Graph_Mapping > Surject_to_bam


#Step2: Surject onto chromosomal paths
rule Surject_to_bam:
    input:
        gam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam",
        gbz = wd + "/temp/{sample}/pangenome.gbz",
    output:
        bam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.bam",
        bam_sort = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam",
        bam_sort_bai = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam.bai",
    log: wd + "/logs/02.surject_to_bam/{sample}."+ref_name+".log"
    benchmark: wd + "/benchmarks/02.surject_to_bam/{sample}."+ref_name+".txt"
    threads: threads
    resources: tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*100
    shell:
        """
        export TMPDIR={wd}/temp/{wildcards.sample}/tmp
        (vg surject --prune-low-cplx -N {wildcards.sample} --read-group 'SM:{wildcards.sample}' -F {wd}/paths.list -x {input.gbz} --interleaved --max-frag-len 3000 -t {threads} -b {input.gam} > {output.bam}) >{log} 2>&1
        samtools view -h {output.bam} -@ {threads} | sed 's/{graph_prefix}chr/chr/g' | samtools sort -@ {threads} - -o {output.bam_sort}
        samtools index -@ {threads} {output.bam_sort}
        # remove the temp pangenome
        rm -rf {wd}/temp/{wildcards.sample}/pangenome*
        """


#Step3.1: Split the bam file by chromosome
rule Bam_preprocess:
    input:
        bam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam",
        chr_size = size_file,
    output:
        bam_split = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.bam",
        bam_RD = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.dedup.bam",
        markdup = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.markdup.txt",
        bam_lefted = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.lefted.bam",
        intervals = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.intervals",
        intervals_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.intervals.bed",
        ex_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.ex.bed"
    log: wd + "/logs/03.bam_processing/{sample}.{chr}."+ref_name+".preprocess.log"
    benchmark: wd + "/benchmarks/03.bam_processing/{sample}.{chr}."+ref_name+".preprocess.txt"
    threads: 6
    resources: tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*5
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        export TMPDIR={wd}
        # Split the bam file by chromosomes and add the RG header line
        samtools view -@ {threads} -h {input.bam} {params.chr} -b -o {output.bam_split}
        samtools index -@ {threads} {output.bam_split}
        # drop duplicate
        (time gatk MarkDuplicates --java-options " -Xms2G -Xmx2G -XX:ParallelGCThreads=6" I={output.bam_split} O={output.bam_RD} M={output.markdup} REMOVE_DUPLICATES=false CREATE_INDEX=true) >> {log} 2>&1
        # Left-alignment
        (cat {output.bam_RD} | bamleftalign -f {ref} > {output.bam_lefted}) >>{log} 2>&1
        # Mark Indel regions for realignment
        samtools index -@ {threads} {output.bam_lefted}
        (/usr/bin/time -v gatk3 -Xms2G -Xmx2G -XX:ParallelGCThreads=6 -T RealignerTargetCreator --remove_program_records --disable_bam_indexing -L {wildcards.chr} -R {ref} -I {output.bam_lefted} -o {output.intervals}) >> {log} 2>&1
        awk -F '[:-]' '{{ if( $3 == "") {{ print ($1"\\t"$2-1"\\t"$2) }} else {{ print ($1"\\t"$2-1"\\t"$3) }} }}' {output.intervals} > {output.intervals_bed}
        bedtools slop -i {output.intervals_bed} -g {input.chr_size} -b 160 > {output.ex_bed}
        """

#Step3.2: Split the bam file by chromosome and realignment
rule Realignment:
    input:
        bam_lefted = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.lefted.bam",
        ex_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.ex.bed",
        ref = ref
    output:
        bam_realign = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.realigned.bam",
    log: wd + "/logs/03.bam_processing/{sample}.{chr}."+ref_name+".realignment.log"
    benchmark: wd + "/benchmarks/03.bam_processing/{sample}.{chr}."+ref_name+".realignment.log"
    threads: 6   
    resources: tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*10
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        # Executing realignment
        rm -rf {wd}/temp/{wildcards.sample}/abra/{params.chr}
        mkdir -p {wd}/temp/{wildcards.sample}/abra/{params.chr}
        (/usr/bin/time -v abra2 --in {input.bam_lefted} --out {output.bam_realign} --ref {ref} --threads {threads} --targets {input.ex_bed} --tmpdir {wd}/temp/{wildcards.sample}/abra/{params.chr}) >> {log} 2>&1
        samtools index -@ {threads} {output.bam_realign}
        """


#Step4: Calling small variants with Deepvariant
rule Deepvariant:
    input:
        bam_realign = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.realigned.bam",
        ref_chr = wd + "/temp/ref/"+ref_name+".{chr}.fa",
    output:
        vcf = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.{chr}.vcf.gz",
        gvcf = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.{chr}.g.vcf.gz",
    log: wd + "/logs/04.deepvariant/{sample}.{chr}."+ref_name+".log"
    benchmark: wd + "/benchmarks/04.deepvariant/{sample}.{chr}."+ref_name+".txt"
    threads: 8
    resources: 
        tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*50,
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        "(/usr/bin/time -v {run_deepvariant} --model_type WGS --ref {input.ref_chr} --reads {input.bam_realign} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --num_shards={threads} --make_examples_extra_args --keep_legacy_allele_counter_behavior=true,--min_mapping_quality=1 ) >{log} 2>&1"

#Step5.1
rule Bam2fastq:
    input:
        bam_split = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.bam",
    output:
        fq_chr = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.fq"
    log: wd + "/logs/05.bam2fastq/{sample}.{chr}."+ref_name+".log"
    benchmark: wd + "/benchmarks/05.bam2fastq/{sample}.{chr}."+ref_name+".txt"
    threads: 4
    resources: 
        tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*5,
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        samtools collate -u -O -@ {threads} {input.bam_split} | samtools fastq -@ {threads} -0 /dev/null > {output.fq_chr}
        """

#Step5.2: PanGenie
rule Pangenie_genotyping:
    input:
        fq_chr = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.fq",
        ref_chr = wd + "/temp/ref/"+ref_name+".{chr}.fa",
        ref_vcf_chr = wd + "/temp/ref/"+ref_name+".{chr}.pangenie.vcf",
    output:
        pangenie_vcf = wd + "/temp/{sample}/pangenie_{sample}."+ref_name+"/{sample}.{chr}_genotyping.vcf.gz",
        pangenie_vcf_index = wd + "/temp/{sample}/pangenie_{sample}."+ref_name+"/{sample}.{chr}_genotyping.vcf.gz.tbi",    
    log: wd + "/logs/05.PanGenie/{sample}.{chr}."+ref_name+".log"
    benchmark: wd + "/benchmarks/05.PanGenie/{sample}.{chr}."+ref_name+".txt"
    threads: 8
    resources: tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=1024*20
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        mkdir -p {wd}/temp/{wildcards.sample}/pangenie_{wildcards.sample}.{ref_name}/{wildcards.sample}
        ({run_pangenie} -i {input.fq_chr} -v {input.ref_vcf_chr} -r {input.ref_chr} -o {wd}/temp/{wildcards.sample}/pangenie_{wildcards.sample}.{ref_name}/{wildcards.sample}.{wildcards.chr} -s {wildcards.sample} -j {threads} -t {threads} -g -e 300000000) > {log} 2>&1
        bgzip -f {wd}/temp/{wildcards.sample}/pangenie_{wildcards.sample}.{ref_name}/{wildcards.sample}.{wildcards.chr}_genotyping.vcf && tabix -f {wd}/temp/{wildcards.sample}/pangenie_{wildcards.sample}.{ref_name}/{wildcards.sample}.{wildcards.chr}_genotyping.vcf.gz
        """

#Step6.1 Joint_calling_dv
rule Glnexus_joint_calling:
    input:
        expand(wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.{chr}.g.vcf.gz",sample=samples,allow_missing=True),
    output:
        small_joint_vcf = wd + "/result/Giraffe_deepvariant.GLnexus.{chr}.vcf.gz"
    log: wd + "/logs/06.glnexus_joint_calling/{chr}."+ref_name+".log"
    benchmark: wd + "/benchmarks/06.glnexus_joint_calling/{chr}."+ref_name+".txt"
    threads: 16 
    resources: tmpdir=wd+"/temp/tmp", mem_mb=1024*20
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        ({run_glnexus} --dir {wd}/temp/tmp_{wildcards.chr} \
            --config DeepVariant \
            {input} \
            --threads {threads} > {wd}/temp/Giraffe_deepvariant.GLnexus.{wildcards.chr}.bcf
            bcftools view {wd}/temp/Giraffe_deepvariant.GLnexus.{wildcards.chr}.bcf -Oz -o {output.small_joint_vcf} && tabix -f {output.small_joint_vcf}) > {log} 2>&1
        """

rule Glnexus_concat_chr:
    input:
        expand(wd + "/result/Giraffe_deepvariant.GLnexus.{chr}.vcf.gz",chr=CHROMOSOMES),
    output:
        wd + "/result/Giraffe_deepvariant.GLnexus.all.vcf.gz"
    log: wd + "/logs/06.glnexus_joint_calling.all."+ref_name+".log"
    benchmark: wd + "/benchmarks/06.glnexus_joint_calling.all."+ref_name+".txt"
    threads: 16 
    resources: tmpdir=wd+"/temp/tmp", mem_mb=1024*20
    shell:
        """
        bcftools concat --threads {threads} $(for i in {1..22} X Y M;do echo {wd}/temp/Giraffe_deepvariant.GLnexus.chr$i.vcf.gz;done) -Oz -o {output}
        tabix -f {output} 
        """

#Step 7.1 pangenie_merge
rule PanGenie_merge:
    input:
        expand(wd + "/temp/{sample}/pangenie_{sample}."+ref_name+"/{sample}.{chr}_genotyping.vcf.gz",sample=samples,allow_missing=True),
    output:
        PanGenie_merge_vcf = wd + "/result/PanGenie.bcftools.{chr}.vcf.gz"
    log: wd + "/logs/07.pangenie_join/{chr}."+ref_name+".log"
    benchmark: wd + "/benchmarks/07.pangenie_join/{chr}."+ref_name+".txt"
    threads: 16 
    resources: tmpdir=wd+"/temp/tmp", mem_mb=1024*20
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        bcftools merge --threads {threads} {input} -Oz -o {output.PanGenie_merge_vcf}
        tabix -f {output.PanGenie_merge_vcf}
        """

rule PanGenie_concat_chr:
    input:
        expand(wd + "/result/PanGenie.bcftools.{chr}.vcf.gz",chr=CHROMOSOMES),
    output:
        wd + "/result/PanGenie.bcftools.all.vcf.gz"
    log: wd + "/logs/07.pangenie_concat.all."+ref_name+".log"
    benchmark: wd + "/benchmarks/07.pangenie_concat.all."+ref_name+".txt"
    threads: 16 
    resources: tmpdir=wd+"/temp/tmp", mem_mb=1024*20
    shell:
        """
        bcftools concat --threads {threads} {input} -Oz -o {output}
        tabix -f {output} 
        """


rule checkends:
    input:
        wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam",
        wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.bam",
        wd + "/result/PanGenie.bcftools.all.vcf.gz",
        wd + "/result/Giraffe_deepvariant.GLnexus.all.vcf.gz"
    output:
        wd + "/logs/{sample}.checkpoints.end"
    threads: 1
    resources: tmpdir= wd + "/temp/tmp/{wildcards.sample}", mem_mb=100
    shell:
        ## check the results
        "touch {output}"
