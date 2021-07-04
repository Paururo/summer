#! /usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import subprocess
import sys
#python3 anotar.py -tsv COV026530.tsv -path /mnt/d/summer/summer/subir/programs/snpEff/
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def parse_args():
    '''
    Parse arguments given to script
    '''
    parser = argparse.ArgumentParser(description = "Parse COVs to get mutations")
    parser.add_argument("-tsv", dest = "tsv_file", required = True)
    parser.add_argument("-path_snpeff", dest = "pathsnpEff", required = True)
    parser.parse_args()
    args = parser.parse_args()
    return args

def clean_VCF(args):
    '''
    Filter VCF file. 
    '''
    import pandas
    vcf = args.tsv_file
    vcf_table = pandas.read_csv(vcf, sep="\t") #read the VCF
    depth_filter = vcf_table['TOTAL_DP'] >= 30   #only positions having more than 30X
    freq_filter  = vcf_table['ALT_FREQ'] >= 0.8 #only alleles bove 80% freq
    vcf_table_filtered = vcf_table[freq_filter & depth_filter].copy() #make a copy of the rows of interest, to avoid warning

    vcf_table_filtered = vcf_table_filtered.sort_values(by=['POS'], ascending = True) #sort data by Position
    vcf_table_filtered.to_csv("filtered_" + vcf, sep="\t", index=False)  

def anotarSnpEff(args):
    vcf = args.tsv_file
    path = args.pathsnpEff
    with open(vcf + '.parannot','w') as out_freqaln:
        with open("filtered_"+vcf,'r') as input_aln:
            contador = 0
            for line in input_aln:
                if contador == 0:
                    linea = '''##fileformat=VCFv4.1\n##contig=<ID=1,length=29903>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'''
                else:
                    line = line.replace("\n","").split("\t")
                    linea = "NC_045512.2\t" + line[1] + "\t.\t" + line[2] + "\t" + line[3] + "\t.\t.\t.\n"
                out_freqaln.write(linea)   
                contador += 1                 
    subprocess.call('java -Xmx4g -jar '+path+'snpEff.jar ann -c ' + path + 'snpEff.config -noStats -no-downstream -no-upstream NC_045512.2 '+ vcf + '.parannot >'+ vcf + '_res.annot', shell = True)

def limpiar(args):
    vcf = args.tsv_file
    with open(vcf + '_clean.annot','w') as out:
        with open(vcf + "_res.annot",'r') as input_vcf:
            titulo = "POS\tREF\tALT\tANNOT\tTYPE\tGENE\tAA\n"
            out.write(titulo)
            for line in input_vcf:
                if "#" not in line:
                    line = line.replace("\n","").split("\t")
                    descripcion = line[7].split("|")
                    linea = line[1] + "\t" + line[3] + "\t" + line[4] + "\t" + descripcion[0] + "\t" + descripcion[1] +"\t" + descripcion[3] +"\t" + descripcion[10] + "\n"
                    out.write(linea)
def main():
    try:
        install(pandas)
    except:
        print("El paquente pandas ya esta instalado")
    args = parse_args()
    clean_VCF(args)
    anotarSnpEff(args)
    limpiar(args)   
if __name__ == '__main__':
    main()