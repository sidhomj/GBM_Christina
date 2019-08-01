import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import os
import pickle
import umap
import matplotlib.pyplot as plt
import seaborn as sns

class JWA(object):
    def __init__(self,Name='Analysis'):
        self.Name = Name

        # Create directory for results of analysis
        directory = self.Name + '_Results'
        self.directory_results = directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Create directory for any temporary files
        directory = self.Name
        if not os.path.exists(directory):
            os.makedirs(directory)

    def Load_Data(self,data_file,id_file,gene_file,mnn_file,Load_Prev_Data=False):
        if Load_Prev_Data is False:
            #Load Data
            pandas2ri.activate()
            readRDS = robjects.r['readRDS']
            X = readRDS(data_file)

            #Load MNN Data
            X_mnn = readRDS(mnn_file)

            #Load Genes
            df_genes = pd.read_csv(gene_file)
            genes = np.array(df_genes.iloc[:,1])

            #Load ID
            df_id = pd.read_csv(id_file)
            cell_id = np.array(df_id.iloc[:,1])

            with open(os.path.join(self.Name,'Data.pkl'),'wb') as f:
                pickle.dump([X,X_mnn,genes,cell_id],f,protocol=4)

        else:
            with open(os.path.join(self.Name, 'Data.pkl'), 'rb') as f:
                X, X_mnn,genes, cell_id = pickle.load(f)

        self.X = X
        self.X_mnn = X_mnn
        self.genes = genes
        self.cell_id = cell_id

    def Run_UMAP(self,Load_Prev_Data=False):
        if Load_Prev_Data is False:
            X_2 = umap.UMAP().fit_transform(self.X_mnn)
            with open(os.path.join(self.Name,'umap_Data.pkl'),'wb') as f:
                pickle.dump(X_2,f,protocol=4)
        else:
            with open(os.path.join(self.Name,'umap_Data.pkl'),'rb') as f:
                X_2 = pickle.load(f)

        self.X_2 = X_2

    def Load_Clustering(self,cluster_file,Load_Prev_Data=False):
        if Load_Prev_Data is False:
            pandas2ri.activate()
            readRDS = robjects.r['readRDS']
            C = readRDS(cluster_file)
            C = np.array(C)

            with open(os.path.join(self.Name,'cluster_data.pkl'),'wb') as f:
                pickle.dump(C,f,protocol=4)

        else:
            with open(os.path.join(self.Name,'cluster_data.pkl'),'rb') as f:
                C = pickle.load(f)

        self.C = C

    def Plot(self,type,gene_name=None,s=15,alpha=1.0,samples=None):
        X_2 = self.X_2
        X = self.X
        sample_id = np.array([x[0:6] for x in self.cell_id])
        C = self.C

        if samples is not None:
            idx = np.where(np.isin(sample_id,samples))[0]
            X_2 = X_2[idx,:]
            sample_id = sample_id[idx]
            X = X[:,idx]
            C = C[idx]

        if type == 'By_Gene':
            c = np.log2(X[np.where(self.genes==gene_name)[0][0],:]+1)
            title = gene_name
            cmap = 'jet'
        elif type == 'By_Sample':
            df = pd.DataFrame()
            df['X'] = X_2[:,0]
            df['Y'] = X_2[:,1]
            df['c'] = sample_id
            title = 'Samples'
        elif type == 'By_Cluster':
            c = C
            c = c.astype(np.str)
            df = pd.DataFrame()
            df['X'] = X_2[:, 0]
            df['Y'] = X_2[:, 1]
            df['c'] = c
            df['c'] = 'Cluster ' + df['c']
            title = 'Clusters'

        plt.figure()
        if type == 'By_Gene':
            plt.scatter(X_2[:,0],X_2[:,1],c=c,s=s,cmap=cmap)
        else:
            sns.scatterplot(data=df, x='X', y='Y', hue='c', linewidth=0, alpha=alpha, s=s)
        plt.xlabel('')
        plt.ylabel('')
        plt.xticks([])
        plt.yticks([])
        plt.title(title)


