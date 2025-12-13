from argparse import ArgumentParser
import sys, os
import numpy as np
import pandas as pd

parser = ArgumentParser(description="De Novo Mutations Discovery in the Offspring, Maher ALnajjar, 2025, MATE-GBI-Genetics and Genomics Dep. Dr.Barta's LAB, Gödöllő, HUNGARY")

parser.add_argument('--depth_min','-d_min', type=int, help="Define the MINIMUM Filtering Depth, default=15", default=15)
parser.add_argument('--depth_max','-d_max', type=int, help="Define the MAXIMUM Filtering Depth, default=45", default=45)
parser.add_argument('--pl_value','-pl_min', type=int, help="Define the MINIMUM Phred-scaled genotype likelihoods (PL) , default=450", default=450)
parser.add_argument('--input_file','-i',help="The file is a standard VCF file coming from GATK for example")
parser.add_argument('--child','-c', type=str, help="list of children sample id(s) separated by comma, otherwise all samples in the vcf will be considered")
parser.add_argument('--output_file','-o', help="output file", default='DNMs_Result')
parser.add_argument('--out_dir','-o_dir', help="indicate an output directory", default='./')
parser.add_argument('--parent1_id','-p1', help="One of the parents id")
parser.add_argument('--parent2_id','-p2', help="The other parent id")

args = parser.parse_args()


def discover_DNMs(vcf_file,output_file,threshold1,threshold2,threshold3,parent1_id,parent2_id,children):
    
    FILE = open(vcf_file,"r")
    
    CHROM=[]
    LOCUS=[]
    REFRENCE=[]
    ALTERATION=[]
    CHILD_GT=[]
    PARENT_1_GT=[]
    PARENT_2_GT=[]
    CHILD_ID=[]
    VARIANT=[]


    for line in FILE:
        if line[0:2] == "#C":
            data = line.rsplit()
            #child_idx = data.index(child_id)
            parent1_idx = data.index(parent1_id)
            parent2_idx = data.index(parent2_id)
            FORMAT_idx = data.index("FORMAT")
            if not children:
                print("all samples in the vcf other than the parents are considered as children")
                children = []
                for chil in data[9:len(data)]:
                    if chil not in [parent1_id,parent2_id]:
                        children.append(chil)
            for child in children:
                globals()[child+"_idx"] = data.index(child)
        elif line[0] != "#":
            data = line.rsplit()
            chr_num = data[0]
            locus = data[1]
            ref = data[3]
            alt = data[4]
            FORMAT = data[FORMAT_idx].split(":")
            DP_idx = FORMAT.index("DP")
            PL_idx = FORMAT.index("PL")
            GQ_idx = FORMAT.index("GQ")
            AD_idx = FORMAT.index("AD")
            
            
            GT_parent1 = data[parent1_idx].split(":")[0]
            GT_parent2 = data[parent2_idx].split(":")[0]
            if (GT_parent1 not in ['0/0','0|0']) or (GT_parent2 not in ['0/0','0|0']):
                #print("skipped, because of  HOMOzygous parents")
                continue
            for child in children:
                GT_COUNT = 0     
                GT_child = data[globals()[child+"_idx"]].split(":")[0]
                if GT_child in ['.','./.','.|.']: 
                    break ## no-calls are not allowed
                if GT_child not in ["0/1","1/0",'0|1','1|0']:
                    continue
                siblings = []
                for CHILD_1 in children:
                    if CHILD_1 != child:
                        siblings.append(CHILD_1)
                for sibling in siblings:
                    globals()["AD_"+sibling] = list(data[globals()[sibling+"_idx"]].split(":"))[AD_idx]
                    globals()["GT_"+sibling] = list(data[globals()[sibling+"_idx"]].split(":"))[0]
                   
                    if (globals()["AD_"+sibling] not in  ["./.","."]) and (globals()["GT_"+sibling] not in ["./.",".",'.|.']):
                        globals()[sibling+"_ad1"] = np.absolute(int(globals()["AD_"+sibling].split(",")[0]))
                        globals()[sibling+"_ad2"] = np.absolute(int(globals()["AD_"+sibling].split(",")[1]))
                        globals()["DP_"+sibling] = int(data[globals()[sibling+"_idx"]].split(":")[DP_idx])
                    else:
                        break
                    if  globals()["GT_"+sibling] in ["0/0","0|0"] and globals()[sibling+"_ad2"] in [0,1]:
                        GT_COUNT +=1
                if GT_COUNT != len(siblings):
                    #print(locus,"skipped, variants found on other sibling")
                    break
                
                else:

                    dp_cou = 0
                    try:

                        DP_child = int(data[globals()[child+"_idx"]].split(":")[DP_idx])
                        PL_child = data[globals()[child+"_idx"]].split(":")[PL_idx]
                        DP_parent1 = int(data[globals()[child+"_idx"]].split(":")[DP_idx])
                        AD_parent1 = data[parent1_idx].split(":")[AD_idx]
                        PL_parent1 = data[parent1_idx].split(":")[PL_idx]
                        DP_parent2 = int(data[parent2_idx].split(":")[DP_idx])
                        AD_parent2 = data[parent2_idx].split(":")[AD_idx]
                        PL_parent2 = data[parent2_idx].split(":")[PL_idx]
                        AD_child = data[globals()[child+"_idx"]].split(":")[AD_idx]
                        parent1_ad1 = np.absolute(int(AD_parent1.split(",")[0]))
                        parent1_ad2 = np.absolute(int(AD_parent1.split(",")[1]))
                        parent2_ad1 = np.absolute(int(AD_parent2.split(",")[0]))
                        parent2_ad2 = np.absolute(int(AD_parent2.split(",")[1]))
                        child_ad1 = np.absolute(int(AD_child.split(",")[0]))
                        child_ad2 = np.absolute(int(AD_child.split(",")[1]))
                        PL_parent1 = np.absolute(int(PL_parent1.split(",")[2]))
                        PL_parent2 = np.absolute(int(PL_parent2.split(",")[2]))
                        PL_child = np.absolute(int(PL_child.split(",")[2]))
                        ad_ratio = 0.0
                    except ValueError:
                        continue
                    if child_ad1 != 0:
                        ad_ratio = float(child_ad2/child_ad1)
                    for sibling in siblings:
                        try:
                            if globals()["DP_"+sibling] >= threshold1 and globals()["DP_"+sibling] <= threshold2:
                                dp_cou  +=1
                        except ValueError:
                            break

                    if len(ref) == len(alt) == 1:
                        variant = "SNP"
                    else:
                        variant = "INDEL" 
                    if parent1_ad2 == 0 and parent2_ad2 == 0 and DP_parent1 >= threshold1  and  DP_parent1 <= threshold2 and DP_parent2 >= threshold1 and DP_parent2 <= threshold2:
                        if DP_child >= threshold1 and DP_child <= threshold2:
                           if True:
                               if child_ad1 > threshold1/2 and child_ad2 > threshold1/2 and ad_ratio >= 0.5  and PL_child >= threshold3 and dp_cou == len(siblings):

                                   CHROM.append(chr_num)
                                   LOCUS.append(locus)
                                   REFRENCE.append(ref)
                                   ALTERATION.append(alt)
                                   CHILD_GT.append(GT_child)
                                   PARENT_1_GT.append(GT_parent1)
                                   PARENT_2_GT.append(GT_parent2)
                                   CHILD_ID.append(child)
                                   VARIANT.append(variant)

    DF = pd.DataFrame(list(zip(CHROM,LOCUS,REFRENCE,ALTERATION,VARIANT,CHILD_GT,PARENT_1_GT,PARENT_2_GT,CHILD_ID)),columns=['CHROM','LOCUS','REFRENCE','ALTERATION','VAR_TYPE','CHILD_GT','PARENT_1_GT','PARENT_2_GT','CHILD_ID'])
    
    DF_sorted = DF.sort_values('CHILD_ID')
    DF_sorted.to_csv(output_file,index=False,sep='\t')

    print("finished!")
    


if __name__ == "__main__":
    if len(sys.argv) == 1:
        os.system("python3 Discover_DNMs.py -h")
    else:
        parent1_num = args.parent1_id
        parent2_num = args.parent2_id
        out_dir = args.out_dir
        denovo_output_file = out_dir +"/"+ args.output_file + ".tsv"
        depth_min = args.depth_min
        depth_max = args.depth_max
        pl_value = args.pl_value
        children = args.child.split(",") if args.child else [] ## Separated by comma
        print("children of interest\t", children)
        
        if args.input_file[-3:] == "vcf":
            discover_DNMs(args.input_file,denovo_output_file,depth_min,depth_max,pl_value,parent1_num,parent2_num,children)
         
        else:
            denovo_raw_file = args.input_file
