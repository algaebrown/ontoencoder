Tranform directed acyclic graph into autoencoder to embed single cell data.

# input: 
- [ddot ontology class](https://ddot.readthedocs.io/en/latest/ontology.html): offers a connection to the ndex, which contains human ontology
- [GO .obo file](http://purl.obolibrary.org/obo/go/go-basic.obo) Seems most of the biological ontology exist in this format
- directed acyclic graph provided by GOATTOOLS

# dependency
- pytorch
- GOATTOOLS
- DDOT
- scanpy

# output:
- variational autoencoder


# function
- select important branches according to set of genes assigned. It is a feature selection problem
   - important branches: ontotypes with high standard devitation
- trim from top (from root) to count layer
- trim from bottom

# training data
- scRNA-seq in AnnData Format. 




