import sys


####Threshold changed from 75 to 80 for calling homozygous and heterozygous alleles


x = sys.argv[1] ##### the filtered vcf file#############
y = sys.argv[2] ####### the prefix for the output file###############

a11 = open(x).readlines()
a22 = [base for base in a11 if not base.startswith("#")] #removal of commented lines

##### Function to analyze the DP4 values ##########################
def DP4analyze(x):
    a2 = x.strip().split(";")
    for each in a2:
        if each.startswith("DP4="):
            a4 = each[4:].strip().split(",")
            R = int(a4[0]) + int(a4[1])
            A = int(a4[2]) + int(a4[3])
            total = R + A
#            print total, R, A
            RF = format(((float(R) / total) * 100), ".2f")
            AF = format(((float(A) / total) * 100), ".2f")
    z = [str(R), str(A), str(RF), str(AF)] + a4 +[str(total)]
    return z

############# Intial filtering based on DP, Base quality and Mapping quality ############
common = "Chromosome Position Ref Var BaseQ Depth MapQ".split(" ")
filter1 = []
failed1 = ["\t".join(common)]
for char in a22:
    b1 = char.strip().split("\t")
    dp = [base[3:] for base in b1[7].split(";") if base.startswith("DP=")]
    mq = [base[3:] for base in b1[7].split(";") if base.startswith("MQ=")]
#    print dp, mq
    if float(b1[5]) >= 50 and float(dp[0]) >= 5 and float(mq[0]) >= 30:
        c1 = b1[:2]+b1[3:6]+ dp+mq+[b1[7]]
        filter1.append("\t".join(c1))
    else:
        c2 = b1[:2]+b1[3:6]+ dp+mq+[b1[7]]
        failed1.append("\t".join(c2))

######### filtering and annotating individual varinat calls ######################################
header = "Chromosome Position Ref Var BaseQ Depth MapQ RefCount AltCount Ref% Alt% Vtype Genotype".split(" ")
filt_ann = []
for char in filter1:
    tmp = char.split("\t")
    dp4 = DP4analyze(tmp[7])
    tmp[5] = dp4[-1]
    if float(dp4[3]) >= 80 and int(dp4[-1]) >5 and int(dp4[6]) >0 and int(dp4[7]) >0:
        if tmp[7].startswith("INDEL"):
            tmp1 = tmp[:7]+dp4[:4]+["INDEL", "Homo"]
            filt_ann.append("\t".join(tmp1))
        else:
            tmp2 = tmp[:7]+dp4[:4]+["SNP", "Homo"]
            filt_ann.append("\t".join(tmp2))
    elif float(dp4[3]) < 75 and int(dp4[-1]) >5 and int(dp4[4]) >0 and int(dp4[5]) >0 and int(dp4[6]) >0 and int(dp4[7]) >0:
        if tmp[7].startswith("INDEL"):
            tmp3 = tmp[:7]+dp4[:4]+["INDEL", "Hetero"]
            filt_ann.append("\t".join(tmp3))
        else:
            tmp4 = tmp[:7]+dp4[:4]+["SNP", "Hetero"]
            filt_ann.append("\t".join(tmp4))
    else:
        failed1.append(char)

############## Sorting the vcf calls ##################################
vcfs = [base.split("\t") for base in filt_ann]
sorted_vcfs = sorted(vcfs, key=lambda x: int(x[1]))
new = ["\t".join(base) for base in sorted_vcfs]

########## identifying the split calls #####################
AltCalls = []
fCalls = []
for each in new:
    each1 = each.strip().split("\t")
    if "," in each1[3]:
        AltCalls.append(each)
    else:
        fCalls.append(each)

################ Writing the output files###############################
final_one = ["\t".join(header)] +fCalls
open(y+"_filtered.vcf","w+").writelines("\n".join(final_one))
open(y+"_failed.vcf","w+").writelines("\n".join(failed1))
if len(AltCalls) >0:
    open(y+"_AltCalls","w+").writelines("\n".join(AltCalls))

