# Fine-mapping

## Description

This a pipeline of summary statistics based fine-mapping using three commonly-used tools.

## Contents

- [Getting started](#4)

- [Usage](#1)
- [Input](#2)
- [Output](#3)

## <a name="4"></a>Getting started

Clone this repo:

```shell
git clone https://github.com/Jianhua-Wang/Fine-mapping.git
```

set up conda environment and download fine-mapping tools:

```shell
cd bin
bash 00_set_up.sh
cd ..
```

build reference panel:

```shell
cd ref
python 01_prepare_reference.py
cd ..
```

:exclamation: In [PAINTOR's framework](https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD), they go through every VCF file when computing LD which is very time-consuming. So I processed the VCF files in [1000G FTP](<ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502>) by splitting them into blocks and convert them into genotype matrix. Although the preparation will cost hours, it can speed up at 30 times. If you cannot wait for the reference preparation, you can get the ready data from `115...212:/f/jianhua/jianhua_pipeline/Fine-mapping/ref`. Also, if you are skilled with python, you can figure out what else I have done with the original VCF files.

## <a name="1"></a>Usage

I wrote the whole pipeline in a python script and it quite simple after activating the conda environment.

I also wrote a brief help:

```shell
(base) ➜  conda activate finemap
(finemap) ➜  python fine_map_pipe.py -h
    
usage: Finemap using three tools, output variants in total credible sets

python fine_map_pipe.py -s 120575 SUM_STAT ./

positional arguments:
  input               input summary statistics
  output              directory of output file, default:working directory

optional arguments:
  -h, --help          show this help message and exit
  -p , --pop          population of input data,[EUR,EAS,SAS,AMR,AFR],
                      default=EUR
  -n , --maxcausal    the maximum number of allowed causal SNPs, default=1
  -s , --samplesize   sample size of input data
  -c , --cred         the cutoff of credible set, default=0.95
  -t , --thread       number of threads, default=10
```

So we can finemap the test data in `input` folder which is can be browsed in [CAUSALdb](<http://mulinlab.tmu.edu.cn/causaldb/block.html?d=957&f=GR018>) using:

```shell
python fine_map_pipe.py -s 120575 ./input/GR018.txt.gz output
```

## <a name="2"></a>Input

Just like the test data, the input file should be a txt file or .txt.gz file containing below information:

| Column | Description                                   |
| ------ | --------------------------------------------- |
| CHR    | chromosome, int, 1-22                         |
| BP     | Base-pair, int, hg19                          |
| rsID   | can be NA                                     |
| MAF    | Minor allele frequency, (0, 0.5]              |
| EA     | Effective allele, forward strand, capital     |
| NEA    | Non-effective allele, forward strand, capital |
| BETA   | Effect size                                   |
| SE     | Standard error                                |
| P      | P value, (0, 1)                               |
| Zscore | Equal to BETA/SE                              |

```shell
zless input/GR018.txt.gz| head -n2
CHR	BP	rsID	MAF	EA	NEA	BETA	SE	P	Zscore
1	762320	exm2268640		C	T	-2.1671804559878427	2.2192489884814903	0.3288	-0.9765377689642318
```

## <a name="3"></a>Output

The output is a txt file including summary statistics and posterior probability of variants in all credible sets.

The example output is:

```shell
head -n2 ./output/GR018_total_credible_set.txt
CHR	BP	rsID	MAF	EA	NEA	BETA	SE	P	Zscore	PAINTOR	CAVIARBF	FINEMAP	block_id	label
1	55505647	rs11591147	0.020899999999999967	T	G	-1.4101773015832226	0.2294923644583753	8.033e-10	-6.144767844069152	0.946675	0.9511994	0.950832	34	7
```

I appended five columns in the input file:

| Column   | Description                                                  |
| -------- | ------------------------------------------------------------ |
| PAINTOR  | Posterior probability computing by PAINTOR                   |
| CAVIARBF | Posterior probability computing by CAVIARBF                  |
| FINEMAP  | Posterior probability computing by FINEMAP                   |
| block_id | No. of block in `/ref/blocks.txt`                            |
| label    | Because the credible set defined by the three tools are different, I used a this column to annotate the type of credible set of the variants. And the rule is PAINTOR:1, CAVIARBF: 2, FINEMAP: 4. So label = 5 means this variant is in the credible set defined by PAINTOR and FINEMAP, but not in CAVIARBF's credible set. |

## Caution

- The default maximum number of allowed causal SNPs is 1, I haven't tested more than 1 causal SNP for many times, so I'm not sure whether it can work in all condition.
- We used the relatively independent LD blocks determined by [ldtect](<https://bitbucket.org/nygcresearch/ldetect-data/src/master/>) which only include three population. And considering the preparation of reference panel, I used blocks of EUR population for all condition as it's the most popular population. 

## Reference

- Schaid, D. J., Chen, W., & Larson, N. B. (2018). From genome-wide associations to candidate causal variants by statistical fine-mapping. *Nature Reviews Genetics*, *19*(8), 491–504. https://doi.org/10.1038/s41576-018-0016-z
- <http://www.christianbenner.com/>
- <https://github.com/gkichaev/PAINTOR_V3.0>
- <https://bitbucket.org/Wenan/caviarbf/src/default/>