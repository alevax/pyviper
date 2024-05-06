from optparse import OptionParser
import pandas as pd
import anndata
import pyviper

optParser = OptionParser()
#optParser.add_option("-n", type="int", dest="num")
optParser.add_option('-n1','--network1',action="store",type = 'string',dest = 'network1')
optParser.add_option('-n2','--network1',action="store",type = 'string',dest = 'network2')
optParser.add_option('-i','--input',action="store",type = 'string',dest = 'input matrix')
optParser.add_option("-o", "--output",action="store",type = 'string',dest="output", default="")


options,args = optParser.parse_args()

ges = anndata.read_csv(options.input, dtype="float64").T
network1 = pd.read_table(options.network1)
network2 = pd.read_table(options.network2)

net_1 = pyviper.Interactome('net2', network1)
net_2 = pyviper.Interactome('net1', network2)
net_1.filter_targets(ges.var_names)
net_2.filter_targets(ges.var_names)

adata_narnea = pyviper.viper(gex_data=ges, # gene expression signature
                            interactome=[net_1,net_2], # list of interactomes
                            eset_filter=False,
                            enrichment = "narnea",
                            njobs=1, # 3 cores
                            verbose=True)

adata_area = pyviper.viper(gex_data=ges, # gene expression signature
                            interactome= [net_1,net_2], # list of interactomes
                            enrichment = "area",
                            njobs=1, # 3 cores
                            verbose=True)



adata_area.to_df().to_csv(options.output +'pyviper_area_nes.csv')
adata_narnea.to_df().to_csv(options.output + 'pyviper_narnea_nes.csv')
pes = pd.DataFrame(adata_narnea.layers['pes'], columns=adata_narnea.to_df().columns, index = adata_narnea.to_df().index)
pes.to_csv(options.output + 'pyviper_narnea_pes.csv')

