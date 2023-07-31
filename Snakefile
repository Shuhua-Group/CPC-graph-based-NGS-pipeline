"""
Author: Lsy & Yzy
Affiliation: HumPOG, Fudan University
Date: June 8th 2023
Aim:  Snakemake workflow for paired-end NGS reads mapping & variant-calling with CPC-graph.
#Requirements: vg== 1.43, samtools, gatk, freebayes(conda), bedtools, abra2, deepvariant
Run:  snakemake  -s Snakefile
Latest modification: July 24th 2023
- todo
"""

# Load configuration from config.yaml
configfile: "config.yaml"
filelist = config['filelist']
wd = config['work_dir']
ref = config['ref']
ref_name = config['ref_name']
gg_file = config['gg_file']
dist_file = config['dist_file']
gbwt_file = config['gbwt_file']
min_file = config['min_file']
xg_file = config['xg_file']
Graph_preprocessed_vcf = config['Graph_preprocessed_vcf']
threads = config['threads']
run_deepvariant = config['run_deepvariant']

# Function to read in the filelist
from functions.basic import parse_filelist
# Function to get the chromsomes names from the graph reference
from functions.basic import get_chromosomes



# Process the filelist
try:
    fq_dict, pair_flag = parse_filelist(filelist)
except (FileNotFoundError, ValueError) as e:
    print(e)

samples = fq_dict.keys()

# Get the chromsomes names
CHROMOSOMES = get_chromosomes(xg_file)


rule all:
    input:
        # final checkpoint for every sample
        expand(wd + "/logs/{sample}.checkpoints.end", sample = samples),
        # merged snv vcf.gz from deepvariant for all samples
        wd + "/result/merged.graph.deepvariant.pass.snp.vcf.gz",
        # merged indel vcf.gz from deepvariant for all samples
        wd + "/result/merged.graph.deepvariant.pass.indel.vcf.gz",
        # merged pangenie vcf.gz for all samples
        wd + "/result/merged.graph.pangenie.vcf.gz",
        # merged vgcall vcf.gz for all samples
        wd + "/result/merged.graph.vgcall.vcf.gz"

#Step0: # Get the chromsomes size
rule get_chromosomes_size:
    input:
        ref
    output:
        chr_size = wd + "/temp/chromosomes_size.txt"
    threads: 1
    resources: mem_mb=1024*1
    shell:
        """
        awk '/^>/{{if (name) print name "\\t" len; name=$0; len=0; next}} {{len+=length($0)}} END{{if (name) print name "\\t" len}}' {input} | sed "s/^>//g" > {output.chr_size}
        """


#Step1: Run vg giraffe for graph reads mapping (~150 CPU hours, 5.5 hours with 36 threads) 
rule Graph_Mapping:
    input:
        fq = lambda wildcards: fq_dict[wildcards.sample]
    output:
        gam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam",
    log: wd + "/logs/01.graph_mapping.{sample}."+ref_name+".log"
    threads: threads
    resources: mem_mb=1024*50
    run:
        # vg giraffe
        cmd_template = "(vg giraffe -d {dist_file} -m {min_file} -x {xg_file} -g {gg_file} -H {gbwt_file} -f "
        fq_cmd = " -f ".join(input.fq)
        shell(cmd_template + fq_cmd + ' > {output.gam}) > {log} 2>&1')
        # vg stat
        cmd = "vg stats -a {output.gam} >> {log}"
        shell(cmd)

# ruleorder: Graph_Mapping > Surject_to_bam


#Step2: Surject onto chromosomal paths
rule Surject_to_bam:
    input:
        gam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam"
    output:
        bam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.bam",
        bam_sort = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam",
        bam_sort_bai = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam.bai",
    log: wd + "/logs/02.surject_to_bam.{sample}."+ref_name+".log"
    threads: threads
    resources: mem_mb=1024*10
    shell:
        """
        # (vg surject --prune-low-cplx {input.gam} > {output.bam}) >{log} 2>&1
        (vg surject --prune-low-cplx -x {xg_file} -b {input.gam} > {output.bam}) >{log} 2>&1
        samtools sort -@ {threads} {output.bam} -o {output.bam_sort}
        samtools index {output.bam_sort}
        """


#Step3.1: Split the bam file by chromosome
rule Bam_preprocess:
    input:
        bam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam",
        chr_size = wd + "/temp/chromosomes_size.txt"
    output:
        bam_split = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.bam",
        bam_lefted = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.lefted.bam",
        intervals = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.intervals",
        intervals_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.intervals.bed",
        ex_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.ex.bed"
    log: wd + "/logs/03.bam_processing.{sample}.{chr}."+ref_name+".preprocess.log"
    threads: 4
    resources: mem_mb=1024*5
    params:
        chr = lambda wildcards: wildcards.chr
    shell:
        """
        # Split the bam file by chromosomes and add the RG header line
        samtools view -@ {threads} -h {input.bam} {params.chr} | awk -v OFS="\\t" -v CHR={params.chr} -v sample={wildcards.sample} 'BEGIN {{header=0}} /^@HD/ {{if (!header) {{print $0; printf("@RG\\tID:" CHR "\\tSM:" sample "\\n"); header=1}} next;}} {{printf("%s\\tRG:Z:" CHR "\\n", $0)}}' | samtools view -S -b -o {output.bam_split} -
        # Left-alignment
        (cat {output.bam_split} | bamleftalign -f {ref} > {output.bam_lefted}) >>{log} 2>&1
        # Mark Indel regions for realignment
        samtools index -@ {threads} {output.bam_lefted}
        (/usr/bin/time -v gatk3 -T RealignerTargetCreator -nt 1 -R {ref} -I {output.bam_lefted} -o {output.intervals}) >> {log} 2>&1
        awk -F '[:-]' '{{ if( $3 == "") {{ print ($1"\\t"$2-1"\\t"$2) }} else {{ print ($1"\\t"$2-1"\\t"$3) }} }}' {output.intervals} > {output.intervals_bed}
        bedtools slop -i {output.intervals_bed} -g {input.chr_size} -b 160 > {output.ex_bed}
        """

#Step3.2: Split the bam file by chromosome and realignment
rule Realignment:
    input:
        bam_split = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.bam",
        bam_lefted = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.lefted.bam",
        intervals = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.intervals",
        intervals_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.intervals.bed",
        ex_bed = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.IndelRealigner.ex.bed",
    output:
        bam_realign = wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.lefted.realigned.bam",
    log: wd + "/logs/03.bam_processing.{sample}.{chr}."+ref_name+".realignment.log"
    threads: threads    
    resources: mem_mb=1024*10
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

#Step3.5: Merge
rule Bam_Merge:
    input:
        expand(wd + "/temp/{sample}/{sample}."+ref_name+".chr_split/{chr}.giraffe.sorted.lefted.realigned.bam",sample = samples,chr = CHROMOSOMES),
    output:
        bam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.lefted.realigned.bam"
    threads: threads
    resources: mem_mb=1024*50
    run:
        bam_now = ["{wd}/temp/{wildcards.sample}/{wildcards.sample}.{ref_name}.chr_split/"+chrom+".giraffe.sorted.lefted.realigned.bam" for chrom in CHROMOSOMES]
        cmd = "samtools merge -@ {threads} {output.bam} " + ' '.join(bam_now) + ' && samtools index -@ {threads} {output.bam}'
        shell(cmd)


#Step4: Calling small variants with Deepvariant
rule Deepvariant:
    input:
        bam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.lefted.realigned.bam"
    output:
        vcf = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.vcf.gz",
        gvcf = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.g.vcf.gz",
    log: wd + "/logs/04.deepvariant.{sample}."+ref_name+".log"
    threads: threads
    resources: mem_mb=1024*50
    shell:
        "(/usr/bin/time -v {run_deepvariant} --model_type WGS --ref {ref} --reads {input.bam} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --num_shards={threads} ) >{log} 2>&1"


#Step5: Calling variants with vg_call
rule vg_call:
    input:
        gam = wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam"
    output:
        vcf = wd + "/result/{sample}/{sample}."+ref_name+".vgcall.vcf.gz",
    log: wd + "/logs/05.vgcall.{sample}."+ref_name+".log"
    threads: threads
    resources: mem_mb=1024*50
    shell:
        """
        (/usr/bin/time -v vg pack -t {threads} -x {xg_file} -g {input.gam} -Q 5 -s 5 -o {wd}/temp/{wildcards.sample}.{ref_name}.pack) >{log} 2>&1
        (/usr/bin/time -v vg call -t {threads} {xg_file} -k {wd}/temp/{wildcards.sample}.{ref_name}.pack -s {wildcards.sample} | bgzip -@ {threads} > {output.vcf}) >>{log} 2>&1
        tabix {output.vcf}
        
        """

ruleorder: Graph_Mapping > vg_call

#Step6: Small variant filtering
rule Vcf_filtering:
    input:
        vcf_d = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.vcf.gz",
        vcf_g = wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.g.vcf.gz"
    output:
        snp = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.snp.vcf",
        indel = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.indel.vcf",
    log: wd + "/logs/06.vcf_filtering.{sample}."+ref_name+".log"
    threads: threads
    resources: mem_mb=1024*20
    shell:
        """
        zcat {input.vcf_d} | grep -E '#|PASS' > {wd}/temp/{wildcards.sample}.graph.deepvariant.pass.vcf
        (gatk SelectVariants -V {wd}/temp/{wildcards.sample}.graph.deepvariant.pass.vcf --select-type-to-include SNP -O {output.snp}) > {log} 2>&1
        (gatk SelectVariants -V {wd}/temp/{wildcards.sample}.graph.deepvariant.pass.vcf --select-type-to-include INDEL -O {output.indel}) >> {log} 2>&1
        """

ruleorder: Deepvariant > Vcf_filtering

#Step7: PanGenie
rule run_Pangenie:
    input:
        snp = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.snp.vcf",
        fq = lambda wildcards: fq_dict[wildcards.sample],
        Graph_preprocessed = f"{Graph_preprocessed_vcf}"
    output:
        reads_combined = wd + "/temp/{sample}/{sample}."+ref_name+".fq",
        pangenie_vcf = wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_genotyping.vcf",
        pangenie_histo = wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_histogram.histo",
        pangenie_fasta = wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_path_segments.fasta",
    log: wd + "/logs/07.panGenie.{sample}."+ref_name+".log"
    threads: threads
    resources: mem_mb=1024*120
    shell:
        """
        zcat {input.fq} > {output.reads_combined}
        mkdir -p {wd}/result/{wildcards.sample}/pangenie_{wildcards.sample}.{ref_name}/{wildcards.sample}
        (/usr/bin/time -v PanGenie -i {output.reads_combined} -v {input.Graph_preprocessed} -r {ref} -o {wd}/result/{wildcards.sample}/pangenie_{wildcards.sample}.{ref_name}/{wildcards.sample} -s {wildcards.sample} -j {threads} -t {threads} -g) > {log} 2>&1
        """

ruleorder: Vcf_filtering > run_Pangenie


#Step6.1 Compress the output vcf
rule Vcf_Compress:
    input:
        snv_vcf = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.snp.vcf",
        indel_vcf = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.indel.vcf",
        pangenie_vcf = wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_genotyping.vcf"
    output:
        snv_vcf = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.snp.vcf.gz",
        indel_vcf = wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.indel.vcf.gz",
        pangenie_vcf = wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_genotyping.vcf.gz"
    threads: threads
    resources: mem_mb=1024*50
    shell:
        """
        bgzip -@ {threads} {input.snv_vcf} && tabix {output.snv_vcf}
        bgzip -@ {threads} {input.indel_vcf} && tabix {output.indel_vcf}
        bgzip -@ {threads} {input.pangenie_vcf} && tabix {output.pangenie_vcf}
        """        

rule checkends:
    input:
        wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.gam",
        wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.bam",
        wd + "/temp/{sample}/{sample}."+ref_name+".giraffe.sorted.bam",
        wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.vcf.gz",
        wd + "/result/{sample}/{sample}."+ref_name+".giraffe.deepvariant.g.vcf.gz",
        wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.snp.vcf.gz",
        wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.indel.vcf.gz",
        wd + "/result/{sample}/{sample}."+ref_name+".vgcall.vcf.gz",
        f"{Graph_preprocessed_vcf}",
        wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_genotyping.vcf.gz"
    output:
        wd + "/logs/{sample}.checkpoints.end"
    threads: 1
    resources: mem_mb=100
    shell:
        ## check the results
        "touch {output}"

# Step8: merge the vcfs
rule merge_results:
    input:
        snv_vcfs = expand(wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.snp.vcf.gz", sample=samples),
        indel_vcfs = expand(wd + "/result/{sample}/{sample}."+ref_name+".graph.deepvariant.pass.indel.vcf.gz", sample=samples),
        pangenie_vcfs = expand(wd + "/result/{sample}/pangenie_{sample}."+ref_name+"/{sample}_genotyping.vcf.gz", sample=samples),
        vgcall_vcfs = expand(wd + "/result/{sample}/{sample}."+ref_name+".vgcall.vcf.gz", sample=samples)
    output:
        final_snp_vcf = wd + "/result/merged.graph.deepvariant.pass.snp.vcf.gz",
        final_indel_vcf = wd + "/result/merged.graph.deepvariant.pass.indel.vcf.gz",
        final_pangenie_vcf = wd + "/result/merged.graph.pangenie.vcf.gz",
        final_vgcall_vcf = wd + "/result/merged.graph.vgcall.vcf.gz",
    threads: threads
    resources: mem_mb=1024*50    
    shell:
        """
        ## merge the results
        bcftools merge --threads {threads} {input.snv_vcfs} | bgzip > {output.final_snp_vcf} && tabix {output.final_snp_vcf}
        bcftools merge --threads {threads} {input.indel_vcfs} | bgzip > {output.final_indel_vcf} && tabix {output.final_indel_vcf}
        bcftools merge --threads {threads} {input.pangenie_vcfs} | bgzip > {output.final_pangenie_vcf} && tabix {output.final_pangenie_vcf}
        bcftools merge --threads {threads} {input.vgcall_vcfs} | bgzip > {output.final_vgcall_vcf} && tabix {output.final_vgcall_vcf}
        """



