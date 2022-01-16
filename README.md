# FastRemap: A Tool for Quickly Remapping Reads between Genome Assemblies 

Remapping tools are used to quickly and efficiently remap genomic data (e.g.,
read data sets) that had been previously mapped to one reference, to another
reference. Remapping tools become more and more relevant with the explosion of
available genetic data sets and references, as traditional methods for mapping
(e.g., read mapping) will have difficulty keeping pace with the computational
requirements. The state-of-the-art remapping tool, CrossMap \[1\], is widely
used but has its shortcomings in performance and end-to-end genome analysis,
that we address with FastRemap. From our evaluation, we find that FastRemap
provides up to a 7.19x speedup, uses as low as 61.7% of the peak memory
consumption, and enables end-to-end downstream analysis compared to CrossMap. 

[\[1\] Zhao et al. *CrossMap: A Versatile Tool for Coordinate Conversion between Genome Assemblies* Bioinformatics 2014.](https://academic.oup.com/bioinformatics/article/30/7/1006/234947?login=true) 

## FastRemap:
    - Currently only supports BAM files as input 
    - if not using gcc 10 or higher, can use the following library: https://github.com/tcbrindle/span

## To clone: 
```
git clone --recurse-submodules git@github.com:CMU-SAFARI/FastRemap.git FastRemap 
```

## To compile:
### zlib: 
```
FastRemap/zlib$ ./configure
FastRemap/zlib$ make
```

### FastRemap: 
may need to use '-lstdc++fs' in LDFLAGS depending on compiler / system. 
```
FastRemap$ make 
```

## To run: 
```
./FastRemap bam [chain file] [input bam file] [unmapped file] [out file]
```
test using the small sample files in test_data folder 
- input / output files should be paths relative to the current directory. 
- e.g., 
	./FastRemap bam test_data/ce6ToCe10.over.chain test_data/little.bam test.unmapped test.out

optional arguments
- --append-tags (-a) to append tags in output bam file 
- --mean (-m) to set insert size 
- --stdev (-s) to set insert_size_stdev
- --times (-t) to set insert_size_fold 

## To validate and compare two bam outputs: 
```
python ./validation/compare_outputs.py [input bam file 1] [input bam file 2] 
```

## To replicate results in the paper: 

Install and run [CrossMap](https://github.com/liguowang/CrossMap) and FastRemap on 
the publicly available DNA-seq read sets:
- human NA12878 illumina read data set (ERR194147 and ERR262997)
- C. elegans N2 illumina read data set (SRR3536210)
- yeast S288C illumina read data set (ERR1938683) 

Using the relevant reference genomes (Genome Sequence Files; \*.fa) and chain files (LiftOver
files; \*.over.chain) at the [UCSC Genome Browser
Site](https://hgdownload.soe.ucsc.edu/downloads.html) 
- human (hg16, hg17, hg18, hg19, hg38) 
- C. elegans (ce2, ce4, ce6, ce10, ce11) 
- yeast (sacCer1, sacCer2, sacCer3) 

For the below example, we will demonstrate the evaluation pipeline on ce10 and ce11. 
First, download the following files: 
- [ce10ToCe11.over.chain](ce10ToCe11.over.chain.gz) 
- [SRR3536210.bam](https://www.ncbi.nlm.nih.gov/sra/SRR3536210) 

First using any read mapping tool, map the SRR3536210.bam file to the old reference genome (e.g., ce10). 
- output to ce10_SRR3536210.bam 

Using the linux time command, get and write all runtime and memory stats output
files to subdirectories ./evaluation/CrossMap and ./evaluation/FastRemap
(depending on which tool was used), e.g.,: 
``` 
/usr/bin/time -v -p -o ./evaluation/FastRemap/ce10_ce11.time FastRemap -c bam ce10ToCe11.over.chain ce10_SRR3536210.bam fastremap_ce10_ce11_unmapped.bed fastremap_ce10_ce11
/usr/bin/time -v -p -o ./evaluation/crossmap/ce10_ce11.time CrossMap.py bam ce10ToCe11.over.chain ce10_SRR3536210.bam ce10_ce11_unmapped.bed > crossmap_ce10_ce11.bam 
``` 

Finally run the plotting script under the evaluation subdirectory: 
```
cd evaluation 
python plot_runtime.py 
```


