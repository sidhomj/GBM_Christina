from JWA.JWA import *

JW_obj = JWA('T')

#T-cells
data_file = 'GBM Single Cell Data Share/072519/resid/t.rds'
mnn_file = 'GBM Single Cell Data Share/072519/mnn/t.rds'
cluster_file = 'GBM Single Cell Data Share/072519/idents/tIdents.RDS'
id_file = 'JW_data/t_id.csv'
gene_file = 'JW_data/t_genes.csv'

JW_obj.Load_Data(data_file,id_file,gene_file,mnn_file,Load_Prev_Data=True)
JW_obj.Run_UMAP(Load_Prev_Data=True)
JW_obj.Load_Clustering(cluster_file,Load_Prev_Data=True)
JW_obj.Cluster_Def(top=10,Load_Prev_Data=True)
# By_Gene, By_Sample, By_Cluster
#samples=['GBM006','GMB030']
JW_obj.Plot(type='By_Gene',gene_name='BATF',s=5)
JW_obj.Plot(type='By_Cluster',alpha=1.0)
JW_obj.Plot(type='By_Sample',samples=samples)
