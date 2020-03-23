import torch
from torch.nn import functional as F
import pandas as pd
import numpy as np
from collections import OrderedDict, defaultdict

class OntoEncoder(torch.nn.Module):
    ''' An gene ontology based variational autoencoder'''
    def __init__(self, topo, gene_name_to_id, device = 'cpu'):
        
        '''instantiate nn.Linear model for every term'''
        super(OntoEncoder, self).__init__()
        self.layers = OrderedDict()
        self.topo = topo
        self.gene_name_to_id = gene_name_to_id
        self.term_name_to_id = pd.Series(np.arange(len(topo)), index = topo.keys())
        self.root = list(self.topo.keys())[-1]
        self.device = device
        
        # initialize linear layer for every term and store in a dictionary
        
        # encoder
        for term in self.topo.keys(): 
            if term != self.root:
                # if the term is not root, create linear layer for it
                self.layers[term, 'encode'] = torch.nn.Linear(len(self.topo[term]['GENES'])+len(self.topo[term]['TERMS']), 1)
        # decoder
        self.reverse_topology() # reverse topology for decoder structure
        for item in self.rev_topo.keys():
            self.layers[item, 'decode'] = torch.nn.Linear(len(self.rev_topo[item]), 1)
    
    def reverse_topology(self):
        ''' to reverse topology for building decoder architecture'''
        self.rev_topo = defaultdict(list)
        for t in self.topo:
            for term in self.topo[t]['TERMS']:
                self.rev_topo[term].append(t) # key as child, value as parents
            for gene in self.topo[t]['GENES']:
                self.rev_topo[gene].append(t)
            
        
    def encode_ont_term(self, x, term):
        ''' encode one term '''
        
        # concatenate the input vector for the term
        if len(self.topo[term]['TERMS']) == 0: # if the term only have genes as child
            input_vector = x[:,self.gene_name_to_id[self.topo[term]['GENES']]]
        elif len(self.topo[term]['GENES']) == 0: # if the term only have terms as child
            input_vector = self.encode_hidden[:,self.term_name_to_id[self.topo[term]['TERMS']]]
        else:
            input_vector = torch.cat([x[:,self.gene_name_to_id[self.topo[term]['GENES']]],self.encode_hidden[:,self.term_name_to_id[self.topo[term]['TERMS']]]],1)
        
        # batch normalization only when input vetor is longer than 1
        if input_vector.shape[1] > 1:
            input_vector = torch.nn.BatchNorm1d(num_features = input_vector.shape[1])(input_vector)
        
        # implement linear --> tanh --> store in hidden layer
        lin_y = self.layers[term, 'encode'](input_vector)
        #tanh = torch.nn.Tanh()
        #y = tanh(lin_y)
        y = F.relu(lin_y)
        self.encode_hidden[:,self.term_name_to_id[term]] = torch.squeeze(y)
            
    def encode(self, x):
        '''comply with the topological order, encode till reaching the root term'''
        self.encode_hidden = torch.zeros(x.shape[0],len(self.topo))
        
        # obey topological order
        for term in self.topo.keys(): 
            if term != self.root:
                self.encode_ont_term(x, term)
        
        # return all terms' hidden state connected to root, this is the embedding.
        z = self.encode_hidden[:,self.term_name_to_id[self.topo[self.root]['TERMS']]]
        return z
    
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5*logvar)
        eps = torch.randn_like(std)
        return mu + eps*std
    
    def decode_ont_term(self, term, is_term = True):
        #print(term)
        input_terms = self.rev_topo[term]
        # concatenate the output vector for the term
        from_embed_layer = [term for term in input_terms if term in self.topo[self.root]['TERMS']]
        other_layer = [term for term in input_terms if term not in self.topo[self.root]['TERMS']]
        #print(len(from_embed_layer), len(other_layer))
        #if len(from_embed_layer) == 0:
        #    x = self.decode_hidden[:, self.term_name_to_id[other_layer]]
        #if len(other_layer) == 0:
        #    x = self.encode_hidden[:, self.term_name_to_id[from_embed_layer]]
        #else:
        x = torch.cat([self.decode_hidden[:, self.term_name_to_id[other_layer]], self.encode_hidden[:, self.term_name_to_id[from_embed_layer]]], 1)
        
        #print(x.shape)
        
        
        #x.unsqueeze_(-1)
        lin_y = self.layers[term, 'decode'](x)
        
        if is_term: # if is term, use transform
            #tanh = torch.nn.Tanh()
            y = F.relu(lin_y)
            self.decode_hidden[:,self.term_name_to_id[term]] = torch.squeeze(y)
        else: # if is gene (output layer, no need to transform)
            y = lin_y
            self.gene_output[:,self.gene_name_to_id[term]] = torch.squeeze(y)
            
    def decode(self, x):
        ''' how should the decoder be like??'''
        self.decode_hidden = torch.zeros(x.shape[0],len(self.topo))
        self.gene_output = torch.zeros(x.shape[0], x.shape[1])
        
        # reverse topological order
        for term in list(self.topo.keys())[:-1][::-1]: 
            if term not in self.topo[self.root]: # if there's a neuron for the term
                self.decode_ont_term(term)
        # gene must be terminal, so can be computed at last
        for gene in self.gene_name_to_id.index:
            self.decode_ont_term(gene, is_term = False)
        
        return self.gene_output
        
    def forward(self, x):
        z = self.encode(x)
        y = self.decode(x)
        return y

################################### Supervised Learning ##################################################
class TopoNet:
    def __init__(self, topo, gene_name_to_id, no_y_class=4, device = 'cpu'):
        '''DCell framework cell type prediction
        input: 
            topo: dict, topology generated by topology.topo_reader
            gene_name_to_id: pandas.Series, index = gene name; values = index in input array X
            no_y_class: int, number of cell type
            
        '''
        super(TopoNet, self).__init__()
        self.layers = OrderedDict()
        self.topo = topo
        self.gene_name_to_id = gene_name_to_id
        self.term_name_to_id = pd.Series(np.arange(len(topo)), index = topo.keys())
        self.root = list(self.topo.keys())[-1]
        self.device = device
        
        # for every term, initialize layer;
        for term in self.topo.keys(): 
            if term == self.root: 
                self.layers[term] = torch.nn.Linear(len(self.topo[term]['GENES'])+len(self.topo[term]['TERMS']), no_y_class)
            else:
                self.layers[term] = torch.nn.Linear(len(self.topo[term]['GENES'])+len(self.topo[term]['TERMS']), 1)
    def forward_ont_term(self, x, term):
        '''forward 1 ontology term
        x: torch.tensor, input data, no_example * no_features
        term: string, 'GO:0005729'
        '''
        if len(self.topo[term]['TERMS']) == 0:
            input_vector = x[:,self.gene_name_to_id[self.topo[term]['GENES']]]
        elif len(self.topo[term]['GENES']) == 0:
            input_vector = self.hidden[:,self.term_name_to_id[self.topo[term]['TERMS']]]
        else:
            input_vector = torch.cat([x[:,self.gene_name_to_id[self.topo[term]['GENES']]],self.hidden[:,self.term_name_to_id[self.topo[term]['TERMS']]]],1)
        
        if input_vector.shape[1] > 1:
            input_vector = torch.nn.BatchNorm1d(num_features = input_vector.shape[1])(input_vector)
        
        
        lin_y = self.layers[term](input_vector)
        
        if term == self.root:
            y = torch.nn.functional.softmax(lin_y)
            
            self.y_pred = y
        else:
            tanh = torch.nn.Tanh()
            y = tanh(lin_y)
            
            self.hidden[:,self.term_name_to_id[term]] = torch.squeeze(y)
                
    def forward(self, x):
        '''comply with the topological order, forward propagate all data
        Input:
            x: torch.tensor, input data, no_example * no_features
        '''
        self.hidden = torch.zeros(x.shape[0],len(self.topo))
        for term in self.topo.keys(): # obey topological order
            self.forward_ont_term(x, term)
        return self.y_pred
            