# PanDepth
<b>PanDepth, an ultra-fast and efficient genomic tool for coverage calculation</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/PanDepth/releases/download/v2.19/PanDepth-2.19-Linux-x86_64.tar.gz
tar zvxf PanDepth-2.19-Linux-x86_64.tar.gz
cd PanDepth-2.19-Linux-x86_64
./pandepth -h
```
### (2) Building from source
```
git clone https://github.com/HuiyangYu/PanDepth.git
cd PanDepth
sh install.sh
cd bin
./pandepth -h
```
## 2. Usage
```
Usage: pandepth -i in.bam [-g gene.gff|-b region.bed] -o outPrefix
 Input/Output options:
   -i    <str>     input of bam/cram file
   -o    <str>     prefix of output file
 Target options:
   -g    <str>     input gff/gtf file for gene region
   -f    <str>     gff/gtf feature type to parse, CDS or exon [CDS]
   -b    <str>     input bed file for list of regions
   -w    <int>     windows size in bp
 Filter options:
   -q    <int>     min mapping quality [0]
   -x    <int>     exclude reads with any of the bits in FLAG set [1796]
 Other options:
   -t    <int>     number of threads [3]
   -r    <str>     input reference genome file for cram decode or GC parse
   -c              enable the calculation of GC content within the target region (requires -r)
   -h              show this help [v2.19]
```
## 3. Example
### 3.1 Perform coverage analysis for each chromosome
```
pandepth -i test.bam -o test1
```
The output file, named "test1.chr.stat.gz", follows the following format:
```
#Chr	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	332615375	331690145	16945298412	99.72	50.95
Chr02	177319215	176853774	8368043949	99.74	47.19
Chr03	289790774	288837978	14143309698	99.67	48.81
Chr04	248932513	248384761	11741149343	99.78	47.17
Chr05	254874144	253897405	12548725234	99.62	49.23
Chr06	253233553	252488900	11854653791	99.71	46.81
Chr07	266382521	265595425	12873116925	99.70	48.33
Chr08	174326481	173874883	8437159467	99.74	48.40
Chr09	278410012	277664711	13068321410	99.73	46.94
```
### 3.2 Perform coverage analysis for each gene
```
pandepth -i test.bam -g test.gff -o test2
```
The output file, named 'test2.gene.stat.gz', follows the following format:
```
#Chr	Start	End	GeneID	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	16621	30840	Caz01g00010.1	819	819	38297	100.00	46.76
Chr01	66931	67440	Caz01g00030.1	510	510	23804	100.00	46.67
Chr01	131205	142847	Caz01g00040.1	3357	3357	161705	100.00	48.17
Chr01	147579	148601	Caz01g00050.1	1023	1023	59443	100.00	58.11
Chr01	165397	166582	Caz01g00060.1	297	297	13764	100.00	46.34
Chr01	180956	190456	Caz01g00070.1	1602	1602	77279	100.00	48.24
Chr01	193092	194509	Caz01g00080.1	411	411	23014	100.00	56.00
Chr01	230152	236238	Caz01g00090.1	1869	1869	94470	100.00	50.55
Chr01	245657	246295	Caz01g00100.1	510	510	29332	100.00	57.51
```
### 3.3 Perform coverage analysis for specific regions
```
pandepth -i test.bam -b test.bed -o test3
```
The output file, named 'test3.bed.stat.gz', follows the following format:
```
#Chr	Start	End	RegionID	Length	CoveredSite	TotalDepth	Coverage(%)	MeanDepth
Chr01	16621	16662	Chr01_16621_16662	42	42	1935	100.00	46.07
Chr01	16838	16947	Chr01_16838_16947	110	110	5522	100.00	50.20
Chr01	17445	17520	Chr01_17445_17520	76	76	3387	100.00	44.57
Chr01	28409	28474	Chr01_28409_28474	66	66	3901	100.00	59.11
Chr01	29747	29917	Chr01_29747_29917	171	171	5141	100.00	30.06
Chr01	30102	30289	Chr01_30102_30289	188	188	8394	100.00	44.65
Chr01	30403	30532	Chr01_30403_30532	130	130	7786	100.00	59.89
Chr01	30805	30840	Chr01_30805_30840	36	36	2231	100.00	61.97
Chr01	66931	67440	Chr01_66931_67440	510	510	23804	100.00	46.67
```
## 4. Speed
The computation time comparison of seven tools for calculating coverage using different numbers of threads when calculating genome coverage with 150Gb of sequencing reads.
![image](https://github.com/HuiyangYu/PanDepth/assets/41780741/36364a18-5e55-4ab7-9daa-115446c64ccb)
## 5. Memory
The Memory requirements comparison of seven tools for calculating coverage using different numbers of threads when calculating genome coverage with 150Gb of sequencing reads.
![image](https://github.com/HuiyangYu/PanDepth/assets/41780741/c7b365a7-d1e0-4ea0-8f6d-4d964f353c45)
## 6. Accuracy
The statistical results of PanDepth on depth and coverage are completely consistent with samtools depth.
## 7. License
-------

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
