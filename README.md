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

The initial release of FastRemap is described in the following paper:
> Jeremie S. Kim, Can Firtina, Meryem Banu Cavlak, Damla Senol Cali, Can Alkan, Onur Mutlu.
> "[**FastRemap: A Tool for Quickly Remapping Reads between Genome Assemblies**](https://arxiv.org/pdf/2201.06255.pdf)".
> Preprint in _arXiv_, Jan 17, 2022.

[\[1\] Zhao et al. *CrossMap: A Versatile Tool for Coordinate Conversion between Genome Assemblies* Bioinformatics 2014.](https://academic.oup.com/bioinformatics/article/30/7/1006/234947?login=true) 

## FastRemap:
    - Currently only supports BAM files as input 
    - if not using gcc 10 or higher, can use the following library: https://github.com/tcbrindle/span

## To clone: 
```
git clone --recurse-submodules https://github.com/CMU-SAFARI/FastRemap.git FastRemap 
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
./FastRemap -f [file type] -c [chain file] -i [input file] -u [output unmapped file] -o [output file]
```

Required arguments: 
- [file_type]:            bam, sam, or bed file depending on input file
- [chain file]:           chain file (https://genome.ucsc.edu/goldenPath/help/chain.html) describes regions of similarity between references
- [input file]:           file containing elements to be remapped based on chain file
- [output unmapped file]: file containing all the elements that couldnt be remapped from the input file based on the provided chain file
- [output file]:          file containing all the remapped elements from the input file

optional arguments
- --append-tags (-a) to append tags in output bam file 
- --mean (-m) to set insert size 
- --stdev (-s) to set insert_size_stdev
- --times (-t) to set insert_size_fold 

BAM test using the small sample files in test_data folder 
- input / output files should be paths relative to the current directory. 
- e.g., 
	./FastRemap -f bam -c test_data/ce6ToCe10.over.chain -i test_data/little.bam -u test.unmapped -o test.out


## To validate and compare two SAM outputs: 
```
python ./validation/compare_outputs.py [input sam file 1] [input sam file 2] 
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

For the below example, we will demonstrate the evaluation pipeline on sacCer1 and sacCer2. 
First, download the following files: 
- [sacCer1ToSacCer2.over.chain](https://hgdownload.soe.ucsc.edu/goldenPath/sacCer1/liftOver/sacCer1ToSacCer2.over.chain.gz) 
- [sacCer1_ERR1938683.bam](https://zenodo.org/record/5945259#.YfyNwRPMI0o) 
	- Alternatively, download an SRA file of sacCer and map it to the old reference genome (e.g., sacCer1), outputting as sacCer1_ERR1938683.bam

Using the linux time command, get and write all runtime and memory stats output
files to subdirectories ./evaluation/CrossMap and ./evaluation/FastRemap
(depending on which tool was used), e.g.,: 
``` 
/usr/bin/time -v -p -o ./evaluation/FastRemap/sacCer1_sacCer2.time FastRemap -f bam -c sacCer1ToSacCer2.over.chain -i sacCer1_ERR1938683.bam -u fastremap_sacCer1_sacCer2_unmapped.bed -o fastremap_sacCer1_sacCer2.bam
/usr/bin/time -v -p -o ./evaluation/crossmap/sacCer1_sacCer2.time CrossMap.py bam sacCer1ToSacCer2.over.chain sacCer1_ERR1938683.bam > crossmap_sacCer1_sacCer2.bam 
``` 

Finally run the plotting script under the evaluation subdirectory: 
```
cd evaluation 
python plot_runtime.py 
```


## Docker Support 

- [DockerHub](https://hub.docker.com/r/alkanlab/fastremap) 
- Dockerfile in top of the tree 

