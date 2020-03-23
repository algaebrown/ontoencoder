import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

########################################### for supervised DCell model ########################################
# read this dataframe for the name of gene ontology terms
term_entropy = pd.read_pickle('/cellar/users/hsher/ontoPrune/data/term_entropy.pickle')

def hidden_to_pandas(model, y_pred, y_true):
    '''
    join hidden state with predicted and true y label
    model: ontoencoder.TopoNet, trained model to provide hidden state
    y_pred: predicted class (not binarized)
    y_true: real class
    
    return pandas.DataFrame, columns = all hidden GO terms and y_pred, y_true, rows = test sample
    '''
    df = pd.DataFrame(model.hidden.data.numpy(), columns = model.term_name_to_id.index)
    df['y_pred'] = y_pred
    df['y_true'] = y_true
    
    return df

def differential_hidden_state(df, cell_group = ['gluta','GABA'], alpha = 0.05, corr_method = 'fdr_bh'):
    ''' use Welch's t-test to find hidden state that are key to distinguish cell type
    df: pandas.DataFrame produced by `hidden_to_pandas`
    cell_group = list of length 2, cell type that your are interested in seeing differential states
    
    return pandas.DataFrame storing result
    '''
    group1 = df.loc[df['y_pred'] == cell_group[0]]
    group2 = df.loc[df['y_pred'] == cell_group[1]]
    stat = [stats.ttest_ind(group1[go_term],group2[go_term], equal_var = False) for go_term in df.columns[:-2]]
    
    # multiple testing correction
    rejected, p_corrected, alpha_Sidak, alpha_bonf = multipletests(pvals = [s[1] for s in stat], alpha = alpha, method = corr_method, is_sorted = False, returnsorted = False)
    
    # store result in dataframe
    stat_df = pd.DataFrame(columns = ['t', 'p-unadjusted', 'reject'])
    stat_df['t'] = [s[0] for s in stat]
    stat_df['p-unadjusted'] = [s[1] for s in stat]
    stat_df['reject'] = rejected
    stat_df.index = df.columns[:-2]

    stat_df['name'] = term_entropy.loc[stat_df.index, 'name']
    stat_df['t_abs'] = stat_df['t'].abs()

    return(stat_df)