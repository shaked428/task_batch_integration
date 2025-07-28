import anndata as ad
import sys
import numpy as np
import pandas as pd
from scib_metrics import bras

## VIASH START
par = {
  'input_integrated': 'resources_test/.../integrated.h5ad',
  'input_solution': 'resources_test/.../solution.h5ad',
  'output': 'output.h5ad',
  'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
  'output': 'output.h5ad',
}
meta = {
  'name': 'bras'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print('Reading input files', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

print('Compute metrics', flush=True)
score = bras(
    X=adata.obsm['X_emb'],
    labels=adata.obs['cell_type'].to_numpy(),
    batch=adata.obs['batch'].to_numpy()
)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['name'] ],
        'metric_values': [ score ]
    }
)

print("Write output AnnData to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
