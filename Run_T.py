from JWA.JWA import *

JW_obj = JWA('T')

#T-cells
data_file = 'GBM Single Cell Data Share/072519/resid/t.rds'
mnn_file = 'GBM Single Cell Data Share/072519/mnn/t.rds'
cluster_file = 'GBM Single Cell Data Share/072519/idents/tIdents.RDS'
id_file = 'JW_data/t_id.csv'
gene_file = 'JW_data/t_genes.csv'

JW_obj.Load_Data(data_file,id_file,gene_file,mnn_file,Load_Prev_Data=True)
alpha_file = 'GBM Single Cell Data Share/072519/vdj/tra.csv'
beta_file = 'GBM Single Cell Data Share/072519/vdj/trb.csv'
JW_obj.Load_TCR(alpha_file,beta_file,Load_Prev_Data=False)
JW_obj.Load_Clustering(cluster_file,Load_Prev_Data=True)
JW_obj.Run_UMAP(Load_Prev_Data=True)
JW_obj.Plot('By_Clone',clone=JW_obj.Clone_Tab.index[0])

list_of_genes = ['IL2','GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','LAG3','TIGIT','PDCD1','HAVCR2','CTLA4']
JW_obj.Cluster_Def(top=50,Load_Prev_Data=True)
list_of_genes = []
for ii,k in enumerate(JW_obj.Cluster_Def_DF.keys(),0):
    list_of_genes.extend(list(JW_obj.Cluster_Def_DF[k]['Gene']))
    # if ii == 0:
    #     break

JW_obj.HM_Clusters(list_of_genes)
JW_obj.Cluster_Prop()
# By_Gene, By_Sample, By_Cluster
#samples=['GBM006','GMB030']
#JW_obj.Plot(type='By_Gene',gene_name='IL7R',s=5)
JW_obj.Plot(type='By_Cluster',alpha=1.0)
JW_obj.Plot(type='By_Sample',samples=samples)

samples=['GBM006','GMB030']
JW_obj.Plot(type='By_Sample',samples=samples)

all_samples = np.unique(np.array([x[0:6] for x in JW_obj.cell_id]))

JW_obj.Plot(type='By_Sample',samples=all_samples[4])
JW_obj.Plot(type='By_Cluster',samples=all_samples[4:])

JW_obj.Plot(type='By_Sample',samples='GMB030')

for s in all_samples:
    JW_obj.Plot(type='By_Sample', samples=s)







