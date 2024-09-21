import anndata as ad
import sys

## VIASH START

par = {
    'input_dataset': 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad',
    'output': 'output.h5ad'
}

meta = {
    'name': 'foo',
    'config': 'bar',
}

## VIASH END

# add helper scripts to path
sys.path.append(meta["resources_dir"])
from utils import _randomize_graph
from read_anndata_partial import read_anndata


print('Read input_dataset', flush=True)
adata = read_anndata(
    par['input_dataset'],
    obs='obs',
    obsp='obsp',
    uns='uns'
)

print('Randomize graph...', flush=True)
adata = _randomize_graph(adata, neighbors_key="knn")

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['name']
adata.write_h5ad(par['output'], compression='gzip')
