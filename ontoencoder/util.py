import pandas as pd
import scanpy as sc
import numpy as np
import anndata
from scipy import stats
import math

def tasic_data():
    ''' load Tasic et.al, Nature data, return as anndata'''
    adata = anndata.read_h5ad('tasic.h5ad')
    return(adata)

def preprocess_X(adata, included_genes):
    ''' preprocess input matrix for training
    1. select only genes included by that ontology
    2. zero pad genes not detected
    3. the input single cell data is log-normalized already
    4. transform to z-score on gene axis
    
    return X, gene_name, gene_name_to_id
    '''
    
    # detected_genes
    detected_genes = list(included_genes.intersection(set(adata.var.index)))
    X = adata[:,detected_genes].X
    
    # zero_padding for undetected genes
    undetected_genes = list(included_genes-set(adata.var.index))
    padding = np.zeros((adata.X.shape[0], len(undetected_genes)))
    X = np.concatenate((X, padding), axis = 1)
    
    # tranform to z score at gene axis
    X = np.nan_to_num(stats.zscore(X, axis = 0),0)
    
    # gene name (column name)
    gene_name = detected_genes + undetected_genes
    gene_name_to_id = pd.Series(np.arange(len(gene_name)),index = gene_name)
    
    return X, gene_name, gene_name_to_id
class DataLoader:
    ''' evenly set data into training set while keeping an eye on balancing label; y needs to be binarized'''
    def __init__(self, X, y):
        self.X = X
        self.y = y
    def split(self, train_size = 0.9, test_size = 0.1):
        # balanced sampler
        
        class_index = {}
        
        for c in range(self.y.shape[1]):
            class_index[c] = np.where(np.where(self.y == 1)[1]== c)[0]
        
        train_index = set()
        # sample for each class
        for c in class_index.keys():
            
            train = np.random.choice(np.array(class_index[c]), math.floor(len(class_index[c])*train_size), replace = False)
            train_index = train_index.union(set(train))
            
        
        test_index = set(np.arange(self.y.shape[0])) - train_index
        
        train_index = list(train_index)
        test_index = list(test_index)
        
        np.random.shuffle(train_index)
        np.random.shuffle(test_index)
        
        return(self.X[train_index,:], self.y[train_index,:], self.X[test_index,:], self.y[test_index,:])