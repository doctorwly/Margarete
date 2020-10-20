import logging,os
import pandas as pd
import numpy as np
import subprocess
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip
import pysam
import argparse
from multiprocessing import Pool

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s') 

def genDict(dim=3):
        if dim==1:
            return defaultdict(int)
        else:
            return defaultdict(lambda: genDict(dim-1))


def getL1Exp(sample,f_bam_file, exp_file,type):
    """
    return combined exp matrix(df)
    """
    logging.info("Generating {sample} expression matrix...".format(sample=sample))
    samfile = pysam.AlignmentFile(f_bam_file, "rb")
    exp = pd.read_csv(exp_file,sep="\t")
    barcodes = exp.columns

    #nesting defaultdicts in an arbitrary depth
    readDict =  genDict()
    n_reads = 0
    for read in samfile:
        last_tag = read.tags[-1]
        if last_tag[0] == "XT":
            geneID = last_tag[1]
            geneID = geneID.split(".")[0]
            attr = read.query_name.split('_')
            if type=="weizhu":
                cell_barcode = attr[1]
                umi = attr[2]
            if type=="cizhu":
                cell_barcode = attr[0]
                umi = attr[1]
            #filter invalid barcode
            if cell_barcode in barcodes: 
                n_reads += 1
                readDict[cell_barcode][umi][geneID] = 1  
    countDict = genDict(dim=2)
    for barcode in readDict:
        for umi in readDict[barcode]:
            for gene in readDict[barcode][umi]:
                countDict[barcode][gene] += 1
    df = pd.DataFrame(countDict)
    combine = pd.concat([df,exp])
    combine = combine.fillna(0)
    combine = combine.astype(int) 
    combine.to_csv("./exp_matrix/{sample}_exp.gz".format(sample=sample),sep="\t",compression='gzip')


def featurecounts(sample,bam_file):
    cmd = """featureCounts -R BAM -F SAF\
 -a /SGRNJ/Database/script/pipe/L1/mytest/data/L1.saf \
 -o featureCounts/{sample} -T 4  -s 0 \
 {bam_file}""".format(sample=sample,bam_file=bam_file)
    os.system(cmd)
    base = os.path.basename(bam_file)
    f_bam = "featureCounts/{base}.featureCounts.bam".format(base=base)
    return f_bam

def run(sample,bam_file,exp_file,type):
    f_bam_file = featurecounts(sample,bam_file)
    getL1Exp(sample,f_bam_file, exp_file,type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('single cell LINE-1')
    parser.add_argument('--mapfile', help='mapfile, 3 columns, "prefix\\tbam\\texp_matrix"', required=True)
    parser.add_argument('--type', help='cizhu or weizhu', required=True)
    args = parser.parse_args()
    mapfile = args.mapfile
    type = args.type
    df_m = pd.read_csv(mapfile,sep="\t",header=None)

    sample_series = df_m.iloc[:,0]
    bam_series = df_m.iloc[:,1]
    exp_series = df_m.iloc[:,2]
    bam_dic = {}
    exp_dic = {}
    for i in range(len(sample_series)):
        bam_dic[sample_series[i]] = bam_series[i]
        exp_dic[sample_series[i]] = exp_series[i]

    os.system("mkdir featureCounts")
    os.system("mkdir exp_matrix")

    p = Pool(3)
    res_list = []
    for sample in sample_series:
        bam_file = bam_dic[sample]
        exp_file = exp_dic[sample]
        res = p.apply_async(run,args=(sample,bam_file,exp_file,type))
        res_list.append(res)
    p.close()
    p.join()
    for res in res_list:
        print (res.get())


