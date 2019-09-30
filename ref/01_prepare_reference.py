#! /usr/bin/env python3

from subprocess import call
import pandas as pd
import os,shutil
from multiprocessing import Pool

pop_list = ['EUR','EAS','SAS','AMR','AFR']

# download 1000G phase3 VCF files

call('mkdir -p ./ld/vcf',shell=True)

for chromosome in range(1,23):
    call(f'wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -P ./ld/vcf',shell=True)
    call(f'wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi -P ./ld/vcf',shell=True)

call('wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -P ./ld/vcf',shell=True)

# get sample id from super population
for pop in pop_list:
    call(f'grep "{pop}" ./ld/vcf/integrated_call_samples_v3.20130502.ALL.panel| cut -f 1 > ./ld/vcf/{pop}.sample',shell=True)
    call(f'mkdir -p ./ld/vcf/{pop}',shell=True)
    call(f'mkdir -p ./ld/vcf/{pop}_1',shell=True)
    call(f'mkdir -p ./ld/txt/{pop}',shell=True)

# split vcf by population and block
blocks = pd.read_csv('./blocks.txt',sep='\t')

def split_reference(pop,chr_id, start, stop):
    call(f'bcftools view -S ./ld/vcf/{pop}.sample -m2 -M2 -i \'INFO/{pop}_AF>0\' ./ld/vcf/ALL.chr{chr_id}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz {chr_id}:{start}-{stop} -O z -o ./ld/vcf/{pop}_1/{pop}_{chr_id}_{start}_{stop}.vcf.gz',shell=True)
    call(f'bcftools norm --rm-dup both ./ld/vcf/{pop}_1/{pop}_{chr_id}_{start}_{stop}.vcf.gz -O z -o ./ld/vcf/{pop}/{pop}_{chr_id}_{start}_{stop}.vcf.gz',shell=True)
    call(f'tabix -f ./ld/vcf/{pop}/{pop}_{chr_id}_{start}_{stop}.vcf.gz',shell=True)

p = Pool(30)
for pop in pop_list:
    if pop in ['EUR','AMR']:
        block = blocks[blocks['pop']=='EUR']
    elif pop in ['EAS','SAS']:
        block = blocks[blocks['pop']=='ASN']
    else:
        block = blocks[blocks['pop']=='AFR']
    for i in blocks.index:
        chr_id, start, stop = blocks.loc[i].values[:-1]
        p.apply_async(split_reference,(pop, chr_id, start, stop))
p.close()
p.join()


# convert vcf to txt (genotype matrix)
def generate_genotype_matrix(pop,chr_id, start, stop):
    vcf = './ld/vcf/{}/{}_{}_{}_{}.vcf.gz'.format(pop,pop,chr_id, start, stop)
    vcf_df = pd.read_csv(vcf,sep='\t',skiprows=256)

    txt_df = pd.DataFrame()
    txt_df[['CHROM','POS','ID','REF','ALT']] = vcf_df[['#CHROM','POS','ID','REF','ALT']]

    def extract_AF(info):
        for ith in info.split(';'):
            if ith.startswith(f'{pop}_AF'):
                return float(ith.split('=')[-1])
    txt_df[f'{pop}_MAF'] = [extract_AF(info) for info in vcf_df['INFO']]
    txt_df[f'{pop}_MAF'] = 0.5-txt_df[f'{pop}_MAF']
    txt_df[f'{pop}_MAF'] = txt_df[f'{pop}_MAF'].abs()
    txt_df[f'{pop}_MAF'] = 0.5-txt_df[f'{pop}_MAF']

    for sample in vcf_df.columns[9:]:
        txt_df[[f'{sample}_1',f'{sample}_2']] = vcf_df[sample].str.split('|',expand=True)

    txt_df = txt_df[txt_df[f'{pop}_MAF']>0]
    txt_df.to_csv('./ld/txt/{}/{}_{}_{}_{}.txt'.format(pop,pop,chr_id, start, stop),sep='\t',index=False)
    call('bgzip -f ./ld/txt/{}/{}_{}_{}_{}.txt'.format(pop,pop,chr_id, start, stop),shell=True)

p = Pool(30)
for pop in pop_list:
    if pop in ['EUR','AMR']:
        block = blocks[blocks['pop']=='EUR']
    elif pop in ['EAS','SAS']:
        block = blocks[blocks['pop']=='ASN']
    else:
        block = blocks[blocks['pop']=='AFR']
    for i in blocks.index:
        chr_id, start, stop = blocks.loc[i].values[:-1]
        p.apply_async(generate_genotype_matrix,(pop,chr_id, start, stop))
p.close()
p.join()

shutil.rmtree('./ld/vcf/')