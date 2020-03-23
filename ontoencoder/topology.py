from collections import OrderedDict, defaultdict
def topo_reader(file):
    '''
    read topology file from DCell format, convert to ordered dictionary
    file: str, path to topology file
    return OrderedDict, keys = [term]['GENES'] or [term]['TERMS'], in topological order (root at last)
    '''
    topo = OrderedDict()
    with open(file) as f:
        x = f.readlines()
        x = [i.rstrip('\n') for i in x]
        
        for i in range(len(x)):
            if 'ROOT:' in x[i]:
                term = x[i].split(' ')[1]
                topo[term] = {}
                
                gene_list = x[i+1].replace('GENES: ','')
                if len(gene_list) > 0:
                    topo[term]['GENES'] = gene_list.split(' ')
                else:
                    topo[term]['GENES'] = []
                    
                term_list = x[i+2].replace('TERMS: ','')
                if len(term_list) > 0:
                    
                    topo[term]['TERMS'] = term_list.split(' ')
                else:
                    topo[term]['TERMS']= []
                
                
    return topo
def included_genes(topo):
    '''returns genes included in topology'''
    all_gene = set()
    for t in topo.keys():
        if 'GENES' in topo[t].keys():
            all_gene = all_gene.union(set(topo[t]['GENES']))
    return all_gene