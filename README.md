# FastRemap: A Tool for Quickly Remapping Reads between Genome Assemblies 

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
