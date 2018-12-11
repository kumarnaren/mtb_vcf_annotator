import cPickle as pk
import sys

###drastic changes as far as prediction is concerned. Now the program takes in Walker format database and generates the output as previously
###version 7 has improvements to include even the intergenic regions that do not form promoter region for any gene
###Made changes to the script for including double mutations in a codon also 
###The variants are reported in a sorted order and consecutive mutations that cause change in the codon are represented by only one one of them
#############################################################################
file1_refDict = sys.argv[1] #the saved annotated dictionary of the ann ref
with open(file1_refDict) as g1:
    refDict = pk.loads(g1.read())
#print "finished loading the ref dict"

file2_mtdb = sys.argv[2] #the mutation database of snps to compare with
f3 = open(file2_mtdb).readlines()

file7_indels = sys.argv[3] ###mutation database of indels to compare with
f7 = open(file7_indels).readlines()

file3_coord = sys.argv[4] # the gene coordinates file
with open(file3_coord) as r1:
    a1 = r1.readlines()

file4_variants = sys.argv[5] # the variants file with indels appended
f1 = (open(file4_variants).readlines())[1:]

file5_output = sys.argv[6] # the prefix for the outfile

file6_gnames = sys.argv[7] # file containing the RvIds and gene names

f8 = sys.argv[8] #outputfolder

###Usage: script.py refDict database coord variants out_prefix
#############################################################################
codontable = {"ATT":"I", "ATC":"I", "ATA": "I", "CTT": "L", "CTC" : "L", "CTA" : "L", "CTG": "L", "TTA": "L", "TTG" :"L",\
              "GTT":"V", "GTC" : "V", "GTA": "V", "GTG" : "V", "TTT": "F", "TTC":"F", "ATG":"M", "TGT":"C","TGC":"C",\
              "GCT":"A", "GCC": "A", "GCA":"A", "GCG":"A", "GGT": "G", "GGC":"G", "GGA": "G", "GGG": "G", "CCT": "P",\
              "CCC": "P", "CCA": "P", "CCG": "P", "ACT":"T", "ACC": "T", "ACA":"T", "ACG":"T", "TCT":"S", "TCC": "S",\
              "TCA": "S", "TCG": "S", "AGT":"S", "AGC": "S", "TAT": "Y", "TAC":"Y", "TGG": "W", "CAA": "Q", "CAG":"Q",\
              "AAT":"N", "AAC":"N", "CAT":"H", "CAC":"H", "GAA":"E", "GAG":"E", "GAT":"D", "GAC":"D", "AAA":"K", "AAG":"K",\
              "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R", "TAA":"*", "TAG":"*", "TGA":"*"}
############################################################################
dlcord = {}
vcord = {}
for each in f1:
    if "SNP" in each:
        vcord[int(each.strip().split("\t")[1])] = each
    else:
        dlcord[int(each.strip().split("\t")[1])] = each

new_list = []
for i in range(1,len(a1),2):
    b1 = a1[i].strip().split("\t")
    b2 = a1[i+1].strip().split("\t")
    b3 = a1[i-1].strip().split("\t")
    if "-" in b2[-1]:
        new = "\t".join([b1[1], b3[3], b1[3], b1[4]])
        new_list.append(new)
    else:
        new1 = "\t".join([b1[1], b3[3], b2[2], b1[4]])
        new_list.append(new1)
#open("tmp1", "wb").writelines("\n".join(new_list))
#print new_list
gene_coord = dict((base.split("\t")[3], base.split("\t")[1:3]) for base in new_list)
#print len(gene_coord), gene_coord["pks1/"], gene_coord["Rv1482c/"]
def vcord_gene(x):
    varDict2 = {}
    for g in gene_coord:
        g1 = gene_coord[g]
        tmp1 = []
        for each in x:
            if each >= int(g1[0]) and each <= int(g1[1]):
                tmp1.append(each)
        varDict2[g] = tmp1
    return varDict2

indelDict1 = vcord_gene(dlcord)
indelDict = dict((base, indelDict1[base]) for base in indelDict1 if len(indelDict1[base]) > 0)
#print indelDict

varDict1 = vcord_gene(vcord)
varDict = dict((base, sorted(varDict1[base])) for base in varDict1 if len(varDict1[base]) > 0)
#print varDict
##############################################################################################
def aminoext(x, y):
    """x is the changed amino acid and y is the list"""
    z = ""
    for i in range(0, len(y), 3):
        if y[i+1] == x:
            z = y[i+2]
    return z
#x1 = aminoext("T", ["3","A","L","3","C","F","3","T","F"])
#print x1
#################################################################################
def newCodon(x, y, t, z,m):
    """x should be the wild type codon and t, y should be the position for change and z, m should be the change"""
    x1 = list(x)
    if t == "-":
        x1[y-1] = z
        codon = x1
    else:
        x1[y-1] = z
        x1[t-1] = m
        codon = x1
    return "".join(codon)
#print newCodon("CCG", 3 ,2,"C", "A")
##############################################################################################
compDict = {"A":"T", "T":"A", "C":"G","G":"C"}
def annVariants(x,y):
    """x should represent the list of annotations for a particular gene
    y denotes the list of coordinates from the vcf file for a gene"""
    ann = []
    ab1 = {int(base.strip().split("\t")[0]):base.strip().split("\t") for base in x}
#    print ab1
    i = 0
    while i < len(y):
        if y[i] in ab1:
            xy1 = ab1[y[i]]
            b1 = (vcord[y[i]]).strip().split("\t")
            if len(xy1) < 9 and xy1[4] == "+":
                b2 = b1 + [".",".", b1[2] +"/"+ xy1[3] +"/"+ b1[3], xy1[4], xy1[5]]
                ann.append("\t".join(b2))
                i+= 1
            elif len(xy1) < 9 and xy1[4] == "-":
                b3 = b1 + [".",".", xy1[2] +"/"+ xy1[3] +"/"+ compDict[b1[3]], xy1[4], xy1[5]]
                ann.append("\t".join(b3))
                i += 1
            elif len(xy1) > 9 and xy1[17] == "+":
                if i < (len(y)-1) and y[i+1] in ab1:
                    xy2 = ab1[y[i+1]]
                    ab2 = (vcord[y[i+1]]).strip().split("\t")
                    if xy2[4] == xy1[4] and y[i+1] - y[i] < 3:
                        a1 = b1[3]
#                        xy2[1] = b1[1]+"|"+xy2[0]
#                        xy2[2] = b1[2]+"|"+vcord[y[i+1]].split("\t")[2]
#                        xy2[3] = b1[3]+"|"+vcord[y[i+1]].split("\t")[3]
                        c1 = newCodon(xy1[4],int(xy1[5]),int(xy2[5]),a1,ab2[3])
                        b4 = b1 + [xy1[4]+"/"+c1, xy1[5], xy2[7]+ "/"+ xy1[3]+"/"+ codontable[c1]]+xy1[-3:-1]
                        b41 = ab2 + [xy1[4]+"/"+c1, xy1[5], xy2[7]+ "/"+ xy1[3]+"/"+ codontable[c1]]+xy1[-3:-1]
                        ann.append("\t".join(b4))
                        ann.append("\t".join(b41))
#                        y.remove(y[i+1])
                        i += 2
                    else:
                        c1 = newCodon(xy1[4], int(xy1[5]), "-", b1[3], "-")
                        b5 = b1 + [xy1[4] + "/" + c1, xy1[5], xy1[7] + "/" + xy1[3] + "/" + codontable[c1]]+xy1[-3:-1]
                        ann.append("\t".join(b5))
                        i += 1
                else:
                    c1 = newCodon(xy1[4], int(xy1[5]), "-", b1[3], "-")
                    b6 = b1 + [xy1[4] + "/" + c1, xy1[5], xy1[7] + "/" + xy1[3] + "/" + codontable[c1]]+xy1[-3:-1]
                    ann.append("\t".join(b6))
                    i += 1
            elif len(xy1) > 9 and xy1[17] == "-":
                if i < (len(y)-1) and y[i+1] in ab1:
                    xy21 = ab1[y[i+1]]
                    ab21 = (vcord[y[i+1]]).strip().split("\t")
                    if xy21[4] == xy1[4] and y[i+1] - y[i] < 3:
                        a2 = b1[3]
   #                     b1[1] = b1[1] + "|" + xy21[0]
   #                     b1[2] = b1[2] + "|" + vcord[y[i + 1]].split("\t")[2]
   #                     b1[3] = b1[3] + "|" + vcord[y[i + 1]].split("\t")[3]
                        c1 = newCodon(xy1[4],int(xy1[5]),int(xy21[5]),compDict[a2],compDict[ab21[3]])
                        b7 = b1 + [xy1[4]+"/"+c1, xy1[5], xy1[7]+ "/"+ xy1[3]+"/"+ codontable[c1]]+xy1[-3:-1]
                        b71 = ab21 + [xy1[4]+"/"+c1, xy1[5], xy1[7]+ "/"+ xy1[3]+"/"+ codontable[c1]]+xy1[-3:-1]
                        ann.append("\t".join(b7))
                        ann.append("\t".join(b71))
   #                     y.remove(y[i+1])
                        i += 2
                    else:
                        c1 = newCodon(xy1[4], int(xy1[5]), "-", compDict[b1[3]], "-")
                        b8 = b1 + [xy1[4] + "/" + c1, xy1[5], xy1[7] + "/" + xy1[3] + "/" + codontable[c1]]+xy1[-3:-1]
                        ann.append("\t".join(b8))
                        i += 1
                else:
                    c1 = newCodon(xy1[4], int(xy1[5]),"-", compDict[b1[3]],"-")
                    b9 = b1 + [xy1[4] + "/" + c1, xy1[5], xy1[7]+ "/" + xy1[3] + "/" + codontable[c1]]+xy1[-3:-1]
                    ann.append("\t".join(b9))
                    i += 1
            else:
                i += 1
        else:
            i += 1
    return ann

def annIndel(x,y):
    """x should represent the list of annotations for a particular gene
    y denotes the list of coordinates from the vcf file for a gene"""
    ann = []
    for each in x:
        a1 = each.strip().split("\t")
        if int(a1[0]) in y:
            b1 = (dlcord[int(a1[0])]).strip().split("\t")
            if len(a1) < 9:
                b2 = dlcord[int(a1[0])].strip("\n").split("\t") + [".",".", b1[2] +"/"+ a1[3] +"/"+ b1[3], a1[4], a1[5]]
                ann.append("\t".join(b2))
            elif len(a1) >9:
                b4 = dlcord[int(a1[0])].strip("\n").split("\t") + [".",".", b1[2] +"/"+ a1[19] +"/"+ b1[3], a1[17], a1[18]]
                ann.append("\t".join(b4))
    return ann
##############################################################################################
#print " annotating the SNPs"
final = []
for each in varDict.keys():
    final += annVariants(refDict[each], varDict[each])

#print "annotating the indels"
for each in indelDict.keys():
    final+= annIndel(refDict[each], indelDict[each])

zfinal = list(set(final))
#print zfinal
#print zfinal ###take annotation from here all the variants called are annotated which is in form of a big list of mutations
##############################################################################################
###### adding gene names function ############
g_a1 = [base.strip() for base in open(file6_gnames).readlines()]
def geneNames(x):
    y = x.split("\t")
    gene = y[17]
    new_list = ""
    for each in g_a1:
        char = each.split("\t")
        if gene in char:
            y[17] = each
            new_list = "\t".join(y)
        elif "/promoter" in y[17] and gene.replace("promoter","") in char:
            y[17] = each
            new_list = "\t".join(y)
        elif "/intergenic" in y[17]:
            y.insert(17, "-")
            new_list = "\t".join(y)
    return new_list

ids_final = [geneNames(base) for base in zfinal]


db_dict = {}
for l1 in f3:
    l = l1.strip().split("\t")
    db_dict["\t".join(l[:2]) + "/\t" + l[3]] = "\t".join(l)

#print db_dict.keys()

def genID(x, y):
    """x is the mutations column as represented in the vcf file and y is the complete line as in the vcf file"""
    mut = ""
    if y[17] != "-":
        mut = y[17] + "\t" + x
    else:
        mut = y[18].strip() + "\t" + x
    return mut


vcf_dict = {}
i = 0
#print len(ids_final)
for line in ids_final:
    d1 = line.strip().split("\t")
    if d1[11] == "SNP":
        vcf_dict[genID(d1[15], d1)] = line.strip()
    else:
        d4 = d1[15].split("/")
        d5 = d4[1]+"/"+d4[2].replace(d4[0], "")
        vcf_dict[genID(d5, d1)+"\t"+d1[11]] = line.strip()

predicted = []
plist = []
for each in vcf_dict:
    e1 = (each.split("\t"))[1].split("/")
#    print each, e1
    for char in db_dict:
        x1 = char.split("\t")
#        print x1, x1[2][0], x1[2][-1]
        if len(x1[2].split("/")) == 3 and x1[2][-1] != "X":
            if "\t".join(x1[1:]) == each and e1[0] != e1[-1]:
                predicted.append(
                    vcf_dict[each].strip() + "\t" + db_dict[char].split("\t")[0] + "\t"+"\t".join(db_dict[char].split("\t")[-2:]))
                plist.append(each)
        elif len(x1[2].split("/")) == 3 and x1[2][-1] == "X":
            x3 = x1[1]+"\t"+x1[2][:-1]
            if each.startswith(x3) and e1[0] != e1[-1]:
                predicted.append(
                    vcf_dict[each].strip() + "\t" + db_dict[char].split("\t")[0] + "\t"+"\t".join(db_dict[char].split("\t")[-2:]))
                plist.append(each)
        else:
            x2 = x1[2].split("/")
            mch = "\t".join([x1[1],x2[0]])
            if each.startswith(mch) and e1[0] != e1[-1]:
                predicted.append(vcf_dict[each].strip() + "\t" + db_dict[char].split("\t")[0] + "\t"+"\t".join(db_dict[char].split("\t")[-2:]))
                plist.append(each)

################################################################################################################


for each in f7:
    x1 = each.strip().split("\t")
    if x1[2].upper() == "INDEL":
        for char in vcf_dict:
            if (char.upper()).endswith("INDEL") and x1[1]+"/" in char and "-" not in char:
                x2 = each.strip().split("\t")
                z1 = vcf_dict[char]+"\t"+ "\t".join([x2[0],x2[-2],x2[-1]])
                predicted.append(z1)
                plist.append(char)
    else:
        if "+" not in x1[1]:
            y1 = x1[3].upper().split("+")
            y2 = []
            y4 = []
            for t1 in vcf_dict:
                t2 = t1.split("\t")
                if t2[0] == x1[1] + "/":
                    y2.append(t2[1])
                    y4.append(t1)

            y3 = "\t".join(y2)
            if y1[0] in y3 and y1[1] in y3:
                z2 = "\t".join("- - - - - - - - - - - - - - -".split(" ")+[x1[3],"-",x1[1]+"/","-",x1[0],x1[-2],x1[-1]])
                predicted.append(z2)
                plist += y4
        else:
            b1 = x1[3].upper().split("+")
            c1 = x1[1].split("+")
            y2 = []
            y5 = []
            for t1 in vcf_dict:
                t2 = t1.split("\t")
                if t2[0] == c1[1] + "/" or t2[0] == c1[0] + "/":
                    y2.append(t2[1])
                    y5.append(t1)
            y3 = "\t".join(y2)
            if b1[0] in y3 and b1[1] in y3:
                z3 = "\t".join("- - - - - - - - - - - - - - -".split(" ") + [x1[3], "-", x1[1], "-", x1[0], x1[-2], x1[-1]])
                predicted.append(z3)
                plist += y5


################################################################################################################

dbgene = list(set([base.strip().split("\t")[1] for base in f3]))
dbgene1 = list(set([base.strip().split("\t")[1] for base in f7] + dbgene))

def filterduplicates(x):
    flist = []
    for each in x:
        x1 = each.split("\t")
        if x1[17].replace("/","") in dbgene:
            flist.append(each)
    return list(set(flist))

##############################################################################################

head_pred = "Chromosome Position Ref Var BaseQ Depth MapQ RefCount AltCount Ref% Alt% Vtype Genotype Codon CodonPos AAchange Strand Gene RvID Drug Confidence References"
head_others = "Chromosome Position Ref Var BaseQ Depth MapQ	RefCount AltCount Ref% Alt%	Vtype Genotype Codon CodonPos AAchange Strand Gene RvID"
#print predicted
#print plist, vcf_dict.keys()[:100]
oth = list(set(vcf_dict.keys()) - set(plist))
tmp1 = []
for each in oth:
    tmp1.append(vcf_dict[each])
tmp1.sort(key=lambda x: int(x.split("\t")[1]))

if len(predicted) >0:
    pz = filterduplicates(predicted)
    pz.sort(key= lambda   x: int(x.split("\t")[1]))
    open(f8+"/"+file5_output+"_predicted","w+").writelines("\n".join(["\t".join(head_pred.split(" "))] + pz))

open(f8+"/"+file5_output+"_others","w+").writelines("\n".join(["\t".join(head_others.split(" "))] + tmp1))










