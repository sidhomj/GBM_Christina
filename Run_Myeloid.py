from JWA.JWA import *

JW_obj = JWA('Myeloid')

#T-cells
data_file = 'GBM Single Cell Data Share/072519/resid/myeloid.rds'
mnn_file = 'GBM Single Cell Data Share/072519/mnn/myeloid.rds'
cluster_file = 'GBM Single Cell Data Share/072519/idents/myeloidIdents.RDS'
id_file = 'JW_data/myeloid_id.csv'
gene_file = 'JW_data/myeloid_genes.csv'

JW_obj.Load_Data(data_file,id_file,gene_file,mnn_file,Load_Prev_Data=True)
JW_obj.Run_UMAP(Load_Prev_Data=True)
JW_obj.Load_Clustering(cluster_file,Load_Prev_Data=True)
# By_Gene, By_Sample, By_Cluster
samples=['GBM006','GMB030']
#JW_obj.Plot(type='By_Gene',gene_name='CD4',samples=samples)
JW_obj.Plot(type='By_Cluster',alpha=1.0,samples=samples)
JW_obj.Plot(type='By_Sample',samples=samples)
