from JWA.JWA import *
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import seaborn as sns

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

tumor_cell_prop.index = tumor_cell_prop.index.str[-3:]
myeloid_cell_prop.index = myeloid_cell_prop.index.str[-3:]

keep = np.intersect1d(tumor_cell_prop.index,myeloid_cell_prop.index)

tumor_cell_prop = tumor_cell_prop.loc[keep]
myeloid_cell_prop = myeloid_cell_prop.loc[keep]

corr = []
p_val = []
tumor_clusters = []
m_clusters = []
for c_t in tumor_cell_prop.columns:
    for c_m in myeloid_cell_prop.columns:
        corr_temp,p = spearmanr(tumor_cell_prop[c_t],myeloid_cell_prop[c_m])
        corr.append(corr_temp)
        p_val.append(p)
        tumor_clusters.append(c_t)
        m_clusters.append(c_m)

_,p_val_corrected,_,_ = multipletests(p_val,alpha=0.10,method='fdr_bh')

df_out = pd.DataFrame()
df_out['Tumor'] = tumor_clusters
df_out['M'] = m_clusters
df_out['Corr'] = corr
df_out['P_Val'] = p_val
df_out['P_Val_Corrected'] = p_val_corrected
df_out.sort_values(by='Corr',ascending=True,inplace=True)
df_out.to_csv('Pos_Corr.csv',index=False)

df_plot = df_out.iloc[0:5]
for it in df_plot.iterrows():
    plt.figure()
    sns.regplot(tumor_cell_prop[it[1]['Tumor']],myeloid_cell_prop[it[1]['M']])
    plt.xlabel('Tumor_'+str(it[1]['Tumor']))
    plt.ylabel('M_'+str(it[1]['M']))




