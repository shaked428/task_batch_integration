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
from utils import _randomize_features
from read_anndata_partial import read_anndata

print('Read input_dataset', flush=True)
adata = read_anndata(
    par['input_dataset'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)

adata.layers['corrected_counts'] = _randomize_features(adata.X)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['name']
adata.write_h5ad(par['output'], compression='gzip')
