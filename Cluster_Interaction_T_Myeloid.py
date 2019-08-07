from JWA.JWA import *
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
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

#####
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

t_cell_prop.index = t_cell_prop.index.str[-3:]
myeloid_cell_prop.index = myeloid_cell_prop.index.str[-3:]

keep = np.intersect1d(t_cell_prop.index,myeloid_cell_prop.index)

t_cell_prop = t_cell_prop.loc[keep]
myeloid_cell_prop = myeloid_cell_prop.loc[keep]

corr = []
p_val = []
t_clusters = []
m_clusters = []
for c_t in t_cell_prop.columns:
    for c_m in myeloid_cell_prop.columns:
        corr_temp,p = spearmanr(t_cell_prop[c_t],myeloid_cell_prop[c_m])
        corr.append(corr_temp)
        p_val.append(p)
        t_clusters.append(c_t)
        m_clusters.append(c_m)

_,p_val_corrected,_,_ = multipletests(p_val,alpha=0.10,method='fdr_bh')

df_out = pd.DataFrame()
df_out['T'] = t_clusters
df_out['M'] = m_clusters
df_out['Corr'] = corr
df_out['P_Val'] = p_val
df_out['P_Val_Corrected'] = p_val_corrected
df_out.sort_values(by='Corr',ascending=False,inplace=True)
df_out.to_csv('Pos_Corr.csv',index=False)

df_plot = df_out.iloc[0:5]
for it in df_plot.iterrows():
    plt.figure()
    sns.regplot(t_cell_prop[it[1]['T']],myeloid_cell_prop[it[1]['M']])
    plt.xlabel('T_'+str(it[1]['T']))
    plt.ylabel('M_'+str(it[1]['M']))


