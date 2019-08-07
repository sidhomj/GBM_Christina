import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import os
import pickle
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from JWA.utils import *
from scipy.stats import spearmanr, mannwhitneyu
from multiprocessing import Pool, cpu_count
import pydiffmap
import phate
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

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

    def Load_Data(self,data_file,counts_data_file,id_file,gene_file,mnn_file,Load_Prev_Data=False):
        if Load_Prev_Data is False:
            #Load Data
            pandas2ri.activate()
            readRDS = robjects.r['readRDS']
            X = readRDS(data_file)

            #Load Counts
            X_counts = readRDS(counts_data_file)

            #Load MNN Data
            X_mnn = readRDS(mnn_file)

            #Load Genes
            df_genes = pd.read_csv(gene_file)
            genes = np.array(df_genes.iloc[:,1])

            #Load ID
            df_id = pd.read_csv(id_file)
            cell_id = np.array(df_id.iloc[:,1])

            with open(os.path.join(self.Name,'Data.pkl'),'wb') as f:
                pickle.dump([X,X_counts,X_mnn,genes,cell_id],f,protocol=4)

        else:
            with open(os.path.join(self.Name, 'Data.pkl'), 'rb') as f:
                X, X_counts, X_mnn, genes, cell_id = pickle.load(f)

        self.X = X
        self.X_counts = X_counts
        self.X_mnn = X_mnn
        self.genes = genes
        self.cell_id = cell_id

    def Load_TCR(self,alpha_file,beta_file,Load_Prev_Data=False):

        if Load_Prev_Data is False:
            df_a = pd.read_csv(alpha_file)
            c = df_a['barcode'].value_counts()
            keep = np.asarray(c[c == 1].index)
            df_a = df_a[df_a['barcode'].isin(keep)]

            df_b = pd.read_csv(beta_file)
            c = df_b['barcode'].value_counts()
            keep = np.asarray(c[c == 1].index)
            df_b = df_b[df_b['barcode'].isin(keep)]

            df_merge = pd.merge(df_a,df_b,on=['barcode'])

            df_out = pd.DataFrame()
            df_out['barcode'] = df_merge['barcode']
            df_out['Clone_ID'] = df_merge['cdr3_nt_x'] + '_' + df_merge['cdr3_nt_y']

            self.barcode_tcr = np.asarray(df_out['barcode'])
            self.clone_id = np.asarray(df_out['Clone_ID'])
            self.Clone_Tab = df_out['Clone_ID'].value_counts()
            self.Clone_Tab = self.Clone_Tab.to_frame()
            self.Clone_Tab['Sample'] = [self.barcode_tcr[np.where(self.Clone_Tab.index[x]==self.clone_id)[0][0]][:6] for x in range(len(self.Clone_Tab))]

    def Run_UMAP(self,Load_Prev_Data=False,sample=None):
        if Load_Prev_Data is False:
            umap_obj = umap.UMAP()
            if sample is not None:
                idx = np.random.choice(range(len(self.X_mnn)),sample,replace=False)
                umap_obj.fit(self.X_mnn[idx])
            else:
                umap_obj.fit(self.X_mnn)

            X_2 = umap_obj.transform(self.X_mnn)

            with open(os.path.join(self.Name,'umap_Data.pkl'),'wb') as f:
                pickle.dump(X_2,f,protocol=4)
        else:
            with open(os.path.join(self.Name,'umap_Data.pkl'),'rb') as f:
                X_2 = pickle.load(f)

        self.X_2 = X_2

    def Run_DiffMap(self,Load_Prev_Data=False,sample=None,k=64):
        if Load_Prev_Data is False:
            if sample is not None:
                idx = np.random.choice(range(len(self.X_mnn)),sample,replace=False)
                dm = pydiffmap.diffusion_map.DiffusionMap.from_sklearn(n_evecs=2,k=k).fit(self.X_mnn[idx])
            else:
                dm = pydiffmap.diffusion_map.DiffusionMap.from_sklearn(n_evecs=2,k=k).fit(self.X_mnn)

            X_2 = dm.transform(self.X_mnn)

            with open(os.path.join(self.Name,'dm_Data.pkl'),'wb') as f:
                pickle.dump(X_2,f,protocol=4)
        else:
            with open(os.path.join(self.Name, 'dm_Data.pkl'), 'wb') as f:
                X_2 = pickle.load(f)

        self.X_2 = X_2

    def Run_Phate(self,Load_Prev_Data=False,sample=None):
        if Load_Prev_Data is False:
            phate_obj = phate.PHATE()
            if sample is not None:
                idx = np.random.choice(range(len(self.X_mnn)),sample,replace=False)
                phate_obj.fit(self.X_mnn[idx])
            else:
                phate_obj.fit(self.X_mnn)

            X_2 = phate_obj.transform(self.X_mnn)

            with open(os.path.join(self.Name, 'phate_Data.pkl'), 'wb') as f:
                pickle.dump(X_2, f, protocol=4)
        else:
            with open(os.path.join(self.Name, 'phate_Data.pkl'), 'rb') as f:
                X_2 = pickle.load(f)

        self.X_2 = X_2

    def Load_Clustering(self,cluster_file,Load_Prev_Data=False,include_clusters=None,exclude_clusters=None):
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

        if exclude_clusters is not None:
            idx_keep = ~np.isin(C, exclude_clusters)
            C = C[idx_keep]
            self.X = self.X[:,idx_keep]
            self.X_mnn = self.X_mnn[idx_keep]
            self.X_counts = self.X_counts[idx_keep]
            self.cell_id = self.cell_id[idx_keep]

        self.C = C

    def Cluster_Def(self,top=10,type='one_v_all',Load_Prev_Data=False,n_jobs=None):
        if Load_Prev_Data is False:
            if type == 'one_v_all':
                DFs = []
                for c in tqdm(np.unique(self.C)):
                    pos = np.where(self.C==c)[0]
                    neg = np.where(self.C!=c)[0]

                    # if n_jobs is None:
                    #     n_jobs = cpu_count()*2
                    #
                    # p = Pool(n_jobs)
                    # args = list(zip(self.genes,
                    #                  [self.genes]*len(self.genes),
                    #                  [self.X_counts]*len(self.genes),
                    #                  [pos]*len(self.genes),
                    #                  [neg]*len(self.genes)))
                    #
                    # p_val = p.starmap(Gene_Enrichment,args)

                    pos_list = []
                    neg_list = []
                    p_val = []
                    for g in self.genes:
                        g_idx = np.where(self.genes==g)[0][0]
                        pos_val = self.X[g_idx,pos]
                        neg_val = self.X[g_idx,neg]
                        pos_list.append(np.mean(pos_val))
                        neg_list.append(np.mean(neg_val))
                        #_,p = mannwhitneyu(pos_val,neg_val)
                        #p_val.append(p)

                    #_,p_val_corrected,_,_ = multipletests(p_val,alpha=0.10,method='fdr_bh')

                    df = pd.DataFrame()
                    df['Gene'] = self.genes
                    df['Pos'] = pos_list
                    df['Neg'] = neg_list
                    df['FC'] = df['Pos']-df['Neg']
                    #df['P_Val'] = p_val
                    #df['P_Val_Adjusted'] = p_val_corrected
                    df.sort_values(by='FC',ascending=False,inplace=True)
                    #df = df[df['P_Val_Adjusted'] < 0.05]
                    df = df.iloc[0:top]
                    DFs.append(df)
            else:
                DFs = []
                for c in tqdm(np.unique(self.C)):
                    pos = np.where(self.C==c)[0]
                    neg = []
                    for o in np.setdiff1d(np.unique(self.C),c):
                        neg.append(np.where(self.C==o)[0])

                    pos_list = []
                    neg_list = []
                    for g in self.genes:
                        g_idx = np.where(self.genes==g)[0][0]
                        pos_val = self.X[g_idx,pos]
                        pos_list.append(np.mean(pos_val))

                        neg_val = []
                        for n in neg:
                            neg_val.append(np.mean(self.X[g_idx,n]))
                        neg_list.append(np.max(neg_val))

                    df = pd.DataFrame()
                    df['Gene'] = self.genes
                    df['Pos'] = pos_list
                    df['Neg'] = neg_list
                    df['FC'] = df['Pos'] - df['Neg']
                    df.sort_values(by='FC', ascending=False, inplace=True)
                    df = df.iloc[0:top]
                    DFs.append(df)

            out = dict(zip(np.unique(self.C),DFs))

            with open(os.path.join(self.Name,'cluster_def.pkl'),'wb') as f:
                pickle.dump(out,f,protocol=4)

            if not os.path.exists(os.path.join(self.directory_results,'Cluster_Def')):
                os.makedirs(os.path.join(self.directory_results,'Cluster_Def'))

            for c,df in zip(np.unique(self.C),DFs):
                df.to_csv(os.path.join(self.directory_results,'Cluster_Def',str(c)+'.csv'),index=False)
        else:
            with open(os.path.join(self.Name,'cluster_def.pkl'),'rb') as f:
                out = pickle.load(f)

        self.Cluster_Def_DF = out

    def Cluster_Prop(self):
        sample_id = np.array([x[0:6] for x in self.cell_id])
        all_samples = np.unique(sample_id)

        prop = []
        for s in all_samples:
            idx = sample_id==s
            p = []
            for c in  np.unique(self.C):
                p.append(np.sum(self.C[idx] == c) / len(self.C[idx]))
            p = np.asarray(p)
            prop.append(p)

        prop = np.vstack(prop)

        df_out = pd.DataFrame(prop)
        df_out.columns = np.unique(self.C)
        df_out.index = all_samples

        self.Cluster_Prop_DF = df_out
        df_out.to_csv(os.path.join(self.directory_results,'Cluster_Prop.csv'))

    def Corr_Genes(self,list_of_genes,dir='pos'):

        if isinstance(list_of_genes,str):
            list_of_genes = [list_of_genes]

        if not os.path.exists(os.path.join(self.directory_results,'Correlation')):
            os.makedirs(os.path.join(self.directory_results,'Correlation'))

        DFs = []
        for g in list_of_genes:
            idx_g = np.where(self.genes==g)[0][0]
            val_g = self.X[idx_g]
            corr_list = []
            for c in np.setdiff1d(self.genes,g):
                idx_c = np.where(self.genes==c)[0][0]
                val_c = self.X[idx_c]
                corr,_ = spearmanr(val_g,val_c)
                corr_list.append(corr)

            df = pd.DataFrame()
            df['Genes'] = np.setdiff1d(self.genes,g)
            df['Corr'] = corr_list

            if dir == 'pos':
                df.sort_values(by='Corr',ascending=False,inplace=True)
            else:
                df.sort_values(by='Corr',ascending=True,inplace=True)

            df.to_csv(os.path.join(self.directory_results,'Correlation', g+'.csv'),index=False)
            DFs.append(df)

        self.Corr_Genes_DF = dict(zip(list_of_genes,DFs))

    def HM_Clusters(self,list_of_genes):
        X = self.X
        idx = np.where(np.isin(self.genes,list_of_genes))[0]
        X = X[idx,:].T

        # sel = np.random.choice(range(len(X)),10000,replace=False)
        # X = X[sel]
        # sns.clustermap(X,standard_scale=1,cmap='bwr')

        idx_e = []
        for c in np.unique(self.C):
            idx_e.append(np.where(self.C==c)[0])
        idx_e = np.hstack(idx_e)
        X = X[idx_e]
        C = self.C[idx_e]

        df = pd.DataFrame(X)
        df.columns = self.genes[idx]
        color_dict = Generate_Color_Dict(self.C)
        row_colors = [color_dict[x] for x in C]
        sns.clustermap(data=df,cmap='bwr',row_cluster=False,col_cluster=True,
                       row_colors=row_colors,standard_scale=1,z_score=None)
        # for c in np.unique(self.C):
        #     idx_c = np.where(self.C==c)[0]
        #     for x in X:
        #         break
        #
        #     break
        # check=1

    def Plot(self,type,gene_name=None,s=15,alpha=1.0,samples=None,clone=None,title=None):
        X_2 = self.X_2
        plt.figure()
        plt.scatter(X_2[:,0],X_2[:,1])
        ax = plt.gca()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        plt.close()

        X = self.X
        sample_id = np.array([x[0:6] for x in self.cell_id])
        C = self.C
        cell_id = self.cell_id

        if samples is not None:
            idx = np.where(np.isin(sample_id,samples))[0]
            X_2 = X_2[idx,:]
            sample_id = sample_id[idx]
            X = X[:,idx]
            C = C[idx]
            cell_id = cell_id[idx]

        if type == 'By_Gene':
            c = X[np.where(self.genes==gene_name)[0][0],:]
            if title is None:
                title = gene_name
            cmap = 'jet'
        elif type == 'By_Sample':
            df = pd.DataFrame()
            df['X'] = X_2[:,0]
            df['Y'] = X_2[:,1]
            df['c'] = sample_id
            if title is None:
                title = 'Samples'
        elif type == 'By_Cluster':
            c = C
            c = c.astype(np.str)
            df = pd.DataFrame()
            df['X'] = X_2[:, 0]
            df['Y'] = X_2[:, 1]
            df['c'] = c
            df['c'] = 'Cluster ' + df['c']
            if title is None:
                title = 'Clusters'
        elif type == 'By_Clone':
            b = self.barcode_tcr[np.isin(self.clone_id, clone)]
            idx = np.isin(cell_id,b)
            if title is None:
                title = 'Clone'

        plt.figure()
        if type == 'By_Gene':
            plt.scatter(X_2[:,0],X_2[:,1],c=c,s=s,cmap=cmap)
        elif type == 'By_Clone':
            plt.scatter(X_2[~idx,0],X_2[~idx,1],c='grey',s=s,alpha=alpha)
            plt.scatter(X_2[idx,0],X_2[idx,1],c='r',s=s,alpha=alpha)
        else:
            # colors = sns.color_palette('Set3', n_colors=len(np.unique(df['c'])))
            #colors = Generate_Color_Dict(df['c'])
            colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c',
                      '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
                      '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
            if len(np.unique(df['c'])) <= len(colors):
                colors = colors[:len(np.unique(df['c']))]
            else:
                colors = sns.color_palette('Set3', n_colors=len(np.unique(df['c'])))
            sns.scatterplot(data=df, x='X', y='Y', hue='c', linewidth=0, alpha=alpha, s=s,palette=colors)

        plt.xlabel('')
        plt.ylabel('')
        plt.xticks([])
        plt.yticks([])
        plt.title(title)
        plt.xlim(xlim)
        plt.ylim(ylim)




