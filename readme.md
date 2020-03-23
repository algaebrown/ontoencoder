Gene ontology(GO)-based autoencoder for embedding single-cell RNA-seq.

Generally the ontoencoder takes three input: X, y and topology (the gene ontology, or any directed acyclic graph you input)

# How to use this directory
please refer to `notebooks/TopoNet*` for examples of supervised learning; `notebooks/OntoEncoder*` for unsupervised learning.

# How to process my training data
Any single-cell RNA-seq can be log-normalized and saved as .h5ad by `scanpy` package.
The processing step is recorded in `notebooks/GSE71585-single cell.ipynb`

# Is there example data
processed data are stored in `/cellar/users/hsher/ontoencoder/notebooks/tasic.h5ad` (accessible to the Ideker lab)

# How to prepare my ontology
the topology should be stored in the DCell format. [See here for an example](https://github.com/idekerlab/DCell/blob/master/training/Topology/GO:0000790_topology)

This topology file can be converted to OntoEncoder/TopoNet-compatible format using `ontoencoder.topology.topo_reader()`)

# How to trim Gene Ontology(GO) into a smaller tree?
Please refer to OntoPrune [https://github.com/algaebrown/ontoPrune] for more information.

# Can I train TopoNet/OntoEncoder on GPU?
I haven't implement that.

# Dependency
`environment.yml` should be helpful. [Refer to here how to install the same conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
