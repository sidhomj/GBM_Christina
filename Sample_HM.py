from JWA.JWA import *
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import umap
import matplotlib.pyplot as plt
import seaborn as sns

JW_obj = JWA('T')

#T-cells
data_file = 'GBM Single Cell Data Share/072519/resid/t.rds'
mnn_file = 'GBM Single Cell Data Share/072519/mnn/t.rds'
cluster_file = 'GBM Single Cell Data Share/072519/idents/tIdents.RDS'
id_file = 'JW_data/t_id.csv'
gene_file = 'JW_data/t_genes.csv'

JW_obj.Load_Data(data_file,id_file,gene_file,mnn_file,Load_Prev_Data=True)
exclude_clusters = ['14', '15']
JW_obj.Load_Clustering(cluster_file,Load_Prev_Data=True)
JW_obj.Cluster_Prop()

t_cell_prop = JW_obj.Cluster_Prop_DF
df_prop = t_cell_prop

g = sns.clustermap(data=df_prop,standard_scale=1,cmap='bwr')
hm = g.ax_heatmap
hm.set_yticklabels(hm.get_yticklabels(),rotation=0)
plt.savefig('T.png')

JW_obj = JWA('Myeloid')

#Myeloid-cells
data_file = 'GBM Single Cell Data Share/072519/resid/myeloid.rds'
mnn_file = 'GBM Single Cell Data Share/072519/mnn/myeloid.rds'
cluster_file = 'GBM Single Cell Data Share/072519/idents/myeloidIdents.RDS'
id_file = 'JW_data/myeloid_id.csv'
gene_file = 'JW_data/myeloid_genes.csv'

JW_obj.Load_Data(data_file,id_file,gene_file,mnn_file,Load_Prev_Data=True)

JW_obj.Load_Clustering(cluster_file,Load_Prev_Data=False)
JW_obj.Cluster_Prop()

myeloid_cell_prop = JW_obj.Cluster_Prop_DF


df_prop = myeloid_cell_prop

g = sns.clustermap(data=df_prop,standard_scale=1,cmap='bwr')
hm = g.ax_heatmap
hm.set_yticklabels(hm.get_yticklabels(),rotation=0)
plt.savefig('Myeloid.png')

JW_obj = JWA('Tumor')

#Tumor-cells
data_file = 'GBM Single Cell Data Share/072519/resid/tumor.rds'
mnn_file = 'GBM Single Cell Data Share/072519/mnn/tumor.rds'
cluster_file = 'GBM Single Cell Data Share/072519/idents/tumorIdents.RDS'
id_file = 'JW_data/tumor_id.csv'
gene_file = 'JW_data/tumor_genes.csv'

JW_obj.Load_Data(data_file,id_file,gene_file,mnn_file,Load_Prev_Data=True)

JW_obj.Load_Clustering(cluster_file,Load_Prev_Data=False)
JW_obj.Cluster_Prop()

tumor_cell_prop = JW_obj.Cluster_Prop_DF
df_prop = tumor_cell_prop
g = sns.clustermap(data=df_prop,standard_scale=1,cmap='bwr')
hm = g.ax_heatmap
hm.set_yticklabels(hm.get_yticklabels(),rotation=0)
plt.savefig('Tumor.png')


t_cell_prop.index = t_cell_prop.index.str[-3:]
myeloid_cell_prop.index = myeloid_cell_prop.index.str[-3:]
tumor_cell_prop.index = tumor_cell_prop.index.str[-3:]

keep = np.intersect1d(tumor_cell_prop.index,myeloid_cell_prop.index)
keep = np.intersect1d(keep,t_cell_prop.index)

t_cell_prop = t_cell_prop.loc[keep]
t_cell_prop.columns = ['T_'+x for x in t_cell_prop.columns]
tumor_cell_prop = tumor_cell_prop.loc[keep]
tumor_cell_prop.columns = ['Tumor_'+x for x in tumor_cell_prop.columns]

myeloid_cell_prop = myeloid_cell_prop.loc[keep]
myeloid_cell_prop.columns = ['Myeloid_'+x for x in myeloid_cell_prop.columns]


df_all = pd.merge(t_cell_prop,tumor_cell_prop,right_index=True,left_index=True)
df_all = pd.merge(df_all,myeloid_cell_prop,right_index=True,left_index=True)

df_prop = df_all
sns.set(font_scale=1.0)
g = sns.clustermap(data=df_prop,standard_scale=1,cmap='bwr')
hm = g.ax_heatmap
hm.set_yticklabels(hm.get_yticklabels(),rotation=0)
plt.savefig('all.png')










