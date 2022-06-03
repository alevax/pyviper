
print("Running test_Pyther.py...")
# some_file.py
import sys, os
from optparse import OptionParser
# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/Pyther_project/pyther_main/libs/')
path_to_this_file = str(os.path.abspath(os.getcwd()))
sys.path.insert(1, path_to_this_file + "/../../libs/")

from pyther_classes import *
from pyther_fn import *

def getStoufferSig(mat):
    my_sig = mat.sum(axis=0)/sqrt(mat.shape[1]*2)
    return(my_sig)

def getMyDateTimeString():
    my_datetime = datetime.now().strftime("%m%d") +\
        datetime.now().strftime("%Y")[2:4] +\
        datetime.now().strftime("%A")[0:3] +\
        datetime.now().strftime("%H%M%S")
    return(my_datetime)

from datetime import datetime
from math import sqrt


### OPTIONS
# create a OptionParser
# class object
parser = OptionParser()

# add options
parser.add_option("-f", "--file",
                  dest = "filename",
                  help = "write report to FILE",
                  metavar = "FILE")
parser.add_option("-q", "--quiet",
                  action = "store_false",
                  dest = "verbose",
                  default = True,
                  help = "don't print status messages to stdout")
parser.add_option("-g", "--ges",
                  dest = "ges",
                  help = "Gene expression signature")
parser.add_option("-i", "--interactome",
                  dest = "interactome",
                  help = "Interactome TSV")
parser.add_option("-v", "--viper_tsv",
                  dest = "viper",
                  help = "Viper TSV computed in R")
parser.add_option("-n", "--out_name_project",
                  dest = "out_project_name",
                  help = "Output convention for file names: project name")
parser.add_option("-o", "--out_dir",
                  dest = "out_dir",
                  help = "Directory for output files")

(options, args) = parser.parse_args()

print("Loading R VIPER result...")
vip_mat = pd.read_csv(options.viper, sep='\t')

# load interactome
print("Loading interactome...")
# intObj = InteractomefromTSV(dirPath + 'data/stoch-10000_pruned.tsv', 's10k-net')
intObj = InteractomefromTSV(options.interactome, options.out_project_name)

# load ges
print("Loading gene expression signature...")
# gesObj = anndata.read_csv(dirPath + 'data/stoch-10000_ccle-ges.tsv', delimiter = '\t')
gesObj = anndata.read_csv(options.ges, delimiter = '\t')

print("Tranposing GES...")
print("(Pyther assumes rows are samples and columns are genes)")
gesObj = anndata.AnnData.transpose(gesObj)

# run aREA
print("Running aREA...")
nesMat = aREA(gesObj, intObj)

print("Tranposing NES matrix...")
nesMat = pd.DataFrame.transpose(nesMat)
nesMat.columns = gesObj.obs_names

print("Computing correlation between Pyther and VIPER...")
nesMat_sig = getStoufferSig(nesMat)
vip_mat_sig = getStoufferSig(vip_mat)
cor_res = pd.Series.corr(nesMat_sig, vip_mat_sig, method = "spearman")

print("Correlation is " + str(cor_res))

print("Saving results...")
# Opening and Closing a file "MyFile.txt"
# for object name file1.
# file1 = open(options.out_dir + getMyDateTimeString() + "_Pyther_VIPER_corr.txt","a")
file1 = open(options.out_dir + options.out_project_name + "_Pyther_VIPER_corr.txt","a")
file1.writelines(str(cor_res))
file1.close()
# adata_final_filename = os.path.join(options.out_dir, options.out_project_name + "_nesMat.h5ad")
# nesMat.write(adata_final_filename)
# adata_final_filename = os.path.join(options.out_dir, options.out_project_name + "_nesMat.csv")
#nesMat.to_csv(adata_final_filename)

print("Done.")
