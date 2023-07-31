# Graph-based-NGS-pipeline

*Alpha test version*

## Summary
Referring to the paper<sup> [[1]](https://www.nature.com/articles/s41586-023-05896-x)[[2]](https://www.nature.com/articles/s41586-023-06173-7)</sup>  and existing protocol<sup>[[3]](https://github.com/vgteam/vg_wdl)</sup>, we have built a pipeline based on *snakemake* to call/genotype snvs/indels/SVs from multi-sample NGS data in HPC.

## Installation
conda >= 22.9.0 is requeired
```
# clone the repo
git clone https://github.com/Shuhua-Group/CPC-graph-based-NGS-pipeline.git
# create the conda enviroment
conda env create -f environment.yaml
```
**[!] There are many bugs while using conda version of DeepVariant, so that we need to [get DeepVariant separately using Docker or Singularity](https://github.com/google/deepvariant).**

## Quick Start
We can run the pipeline on the toydata for bug test.
```
cd CPC-graph-based-NGS-pipeline
```

### (1). Modify the configfiles
At the root directory of the pipeline, here is `config.yaml` and `toydata/test.filelist`
Replace `{download_dir}` to the absolute directory where the pipeline was downloaded:
`config.yaml`:
```
## filelist
## Two columns for sample names and download_dir of fq files, if paired-end sequenced, the download_dir of two fq files are joint with ';' 
filelist: {download_dir}/toydata/test.filelist

## Working Dir
work_dir: {download_dir}/toydata/out
## Number of threads for a individual at every single step
threads: 4
## Absolute directories to graph reference files
gg_file: {download_dir}/toydata/graph/MHC.gg
...
```
`toydata/test.filelist`
```
HG00403	{download_dir}/toydata/fq/HG00403.R1.fastq.gz;{download_dir}/toydata/fq/HG00403.R2.fastq.gz
NA06984	{download_dir}/toydata/fq/NA06984.R1.fastq.gz;{download_dir}/toydata/fq/NA06984.R2.fastq.gz
NA18486	{download_dir}/toydata/fq/NA18486.R1.fastq.gz;{download_dir}/toydata/fq/NA18486.R2.fastq.gz
NA18525	{download_dir}/toydata/fq/NA18525.R1.fastq.gz;{download_dir}/toydata/fq/NA18525.R2.fastq.gz
```
Add the Absolute directory of executable `run_deepvariant` in DeepVariant to `config.yaml`:
```
## Absolute directory of executable `run_deepvariant` in DeepVariant
run_deepvariant: {dir_to_run_deepvariant}
```


### (2). Run snakemake
Dry-run the pipline at first:
```
snakemake -n
```
Run the pipeline locally, at least 100 GB of memory is required:
```
snakemake -s ./Snakefile --configfile ./config.yaml -c 16
```
You can also run the pipeline in PBS or SLURM system, see details at [snakemake doc](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### (3). Results
If the pipeline runs correctly, the results file will be written to `{download_dir}/toydata/out/result`, including:
```
# merged snv vcf.gz from deepvariant for all samples
merged.graph.deepvariant.pass.snp.vcf.gz

# merged indel vcf.gz from deepvariant for all samples
merged.graph.deepvariant.pass.indel.vcf.gz

# merged pangenie vcf.gz for all samples
merged.graph.pangenie.vcf.gz

# merged vgcall vcf.gz for all samples
merged.graph.vgcall.vcf.gz
```
**Note: the reference in toydata is a partial graph at HLA loci, but it does not mean that the variation sets obtained by running this pipeline on a local graph are correct and accurate**.

## Run on real data
### Pangenome graph:
You can get the whole-genome minigraph-cactus graph from [HPRC](https://github.com/human-pangenomics/hpp_pangenome_resources) and [CPC](https://github.com/Shuhua-Group/Chinese-Pangenome-Consortium-Phase-I/tree/main), and modify the `config.yaml`:
```
## Absolute directories to graph reference files
gg_file: {graph_dir}/CPC.Phase1.CHM13v2-minaf.0.1.gg
dist_file: {graph_dir}/CPC.Phase1.CHM13v2-minaf.0.1.dist
gbwt_file: {graph_dir}/CPC.Phase1.CHM13v2-minaf.0.1.gbwt
min_file: {graph_dir}/CPC.Phase1.CHM13v2-minaf.0.1.min
xg_file: {graph_dir}/CPC.Phase1.CHM13v2-minaf.0.1.xg

## PanGenie reference vcf
Graph_preprocessed_vcf: {graph_dir}/CPC.Phase1.CHM13v2-minaf.0.1.pangenie.vcf
```

### NGS data:
You need to prepare your own filelist: two columns for sample names and absolute directories of fq files, if paired-end sequenced, the directories of two fq files are joint with ';'.  
```
HG00403	{abs_dir}/HG00403.R1.fastq.gz;{abs_dir}/HG00403.R2.fastq.gz
NA06984	{abs_dir}/NA06984.R1.fastq.gz;{abs_dir}/NA06984.R2.fastq.gz
NA18486	{abs_dir}/NA18486.R1.fastq.gz;{abs_dir}/NA18486.R2.fastq.gz
NA18525	{abs_dir}/NA18525.R1.fastq.gz;{abs_dir}/NA18525.R2.fastq.gz
```
In real-data testing, we used a 36-core server to analyse pair-ends ~50x NGS data from three samples, taking a total of ~20 hours and consuming a peak of ~330 GB of memory.
