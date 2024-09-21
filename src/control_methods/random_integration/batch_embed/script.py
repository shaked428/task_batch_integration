import anndata as ad
import sys

## VIASH START

par = {
    'input_dataset': 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad',
    'output': 'output.h5ad',
}

meta = { 
    'functionality': 'foo',
    'config': 'bar'
}

## VIASH END
sys.path.append(meta["resources_dir"])
from utils import _randomize_features
from read_anndata_partial import read_anndata


print('Read input_dataset', flush=True)
adata = read_anndata(
    par['input_dataset'],
    obs='obs',
    obsm='obsm',
    uns='uns'
)

print('Process data...', flush=True)
adata.obsm["X_emb"] = _randomize_features(
    adata.obsm["X_pca"],
    partition=adata.obs["batch"]
)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['name']
adata.write_h5ad(par['output'], compression='gzip')
