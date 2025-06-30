import anndata as ad
from typing import Any, Callable, Optional, Dict
import sys

import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.stats import rankdata
from tqdm import tqdm

def ecdf_qnorm_genes(
    adata: sc.AnnData,
    batch_key: str,
    layer: Optional[str] = None,
    merge_grid_size: int = 200,
    summary_function: Callable = np.nanmean,
    min_num_for_freq: int = 20,
    remove_nan: bool = True,
    min_ref_val: int = 4,
    verbose: bool = True,
) -> sparse.csr_matrix:
    """
    Run quantile normalization each gene, with respect to batches
    :param adata: Input adata
    :param layer: layer name
    :param batch_key: batch key in ``adata.obs``
    :param merge_grid_size: Precision of ecdf, larger -> more accurate but slower
    :param summary_function: function to norm with, such as np.nanmean np.nanmedian, etc
    :param min_num_for_freq:
    :param remove_nan: remove rows with nan values
    :param min_ref_val:
    :param verbose: Whether to print progress
    :return: a csr_matrix with integrated values
    """
    if layer:
        data_in: sparse.spmatrix = adata.layers[layer]
    else:
        data_in: sparse.spmatrix = adata.X

    # Batch initialization
    batch_ids = adata.obs[batch_key].values
    unique_batches = np.unique(batch_ids)
    num_batches = len(unique_batches)
    batch_key_to_index: Dict[Any, np.ndarray] = {  # noqa
        batch: batch_ids == batch for batch in batch_ids
    }

    num_samples, num_genes = data_in.shape
    non_zero_elements = (
        sparse.csr_matrix.count_nonzero(data_in)
        if sparse.issparse(data_in)
        else np.count_nonzero(data_in)
    )
    if num_samples != len(batch_ids):
        raise ValueError("Batch ID length does not match number of samples")

    # Output initialization

    exp_output = np.zeros((num_samples, num_genes))
    merge_grid = np.linspace(0, 1, merge_grid_size + 1)

    if verbose:
        print("Ranking Batch Data:")

    freq_matrix = np.zeros(data_in.shape)
    with tqdm(total=len(batch_key_to_index), disable=not verbose) as pbar:
        for batch_idx, batch in enumerate(list(batch_key_to_index.keys())):
            batch_data = (
                data_in[batch_key_to_index[batch], :].toarray()
                if sparse.issparse(data_in)
                else data_in[batch_key_to_index[batch], :]
            )
            ranked_data = rankdata(batch_data, method="average", axis=0)
            freq_data = ranked_data - np.min(ranked_data, axis=0)
            freq_data = freq_data / np.maximum(np.max(freq_data, axis=0), 1)
            freq_matrix[batch_key_to_index[batch], :] = freq_data
            pbar.update(1)

    data_in = data_in.tocsc() if sparse.issparse(data_in) else data_in

    with tqdm(total=num_genes, disable=not verbose) as pbar:
        for gene_index in range(num_genes):
            cumulative_freq = np.zeros(num_samples)
            merge_grid_matrix = np.full((num_batches, len(merge_grid)), np.nan)

            gene_data = (
                data_in.getcol(gene_index).toarray().flatten()
                if sparse.issparse(data_in)
                else data_in[:, gene_index]
            )
            for batch_index, batch in enumerate(unique_batches):
                batch_data = gene_data[batch_key_to_index[batch]]
                normalized_ranks = freq_matrix[batch_key_to_index[batch], gene_index]
                cumulative_freq[batch_key_to_index[batch]] = normalized_ranks

                unique_batch_data, unique_indices = np.unique(
                    batch_data.data, return_index=True
                )

                if (
                    len(unique_batch_data) > min_num_for_freq
                    and len(unique_indices) >= min_ref_val
                ):
                    x_normalized_ranks = normalized_ranks[unique_indices]
                    x_batch_data = batch_data[unique_indices]

                    merge_grid_matrix[batch_index, :] = np.maximum(
                        np.interp(
                            merge_grid,
                            x_normalized_ranks,
                            x_batch_data,
                            left=np.nan,
                            right=np.nan,
                        ),
                        0,
                    )
            mean_exp = summary_function(merge_grid_matrix, axis=0)
            non_zero_indices = np.nonzero(cumulative_freq)
            exp_output[non_zero_indices, gene_index] = np.interp(
                cumulative_freq[non_zero_indices], merge_grid, mean_exp
            )
            pbar.update(1)

    if remove_nan:
        nan_columns = np.zeros(num_genes, dtype=int)
        for i in range(num_genes):
            nan_mask = np.isnan(exp_output[:, i])
            num_nans = np.sum(nan_mask)
            if num_nans > 0:
                nan_columns[i] = num_nans
                exp_output[nan_mask, i] = 0
        if verbose:
            print(
                f"Total columns with nan values: {np.sum(nan_columns > 0)}, "
                f"Total nan values: {np.sum(nan_columns)}"
            )
    if verbose:
        print(
            f"Done gene quantile matching. Normalized {num_batches} batches.\n"
            f"NNZ before: {non_zero_elements} ("
            f"{non_zero_elements / (num_samples * num_genes)}), "
            f"after: {np.count_nonzero(exp_output)} "
            f"({np.count_nonzero(exp_output) / (num_samples * num_genes)})"
        )

    return sparse.csr_matrix(exp_output)

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
  'output': 'output.h5ad',
  'batch_key': 'batch'
}
meta = {
  'name': 'ecdfQnorm'
}
## VIASH END


print("try only imports, load input and write output")
print('Reading input files', flush=True)
# input = ad.read_h5ad(par['input'])

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('>> Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/counts',
    obs='obs',
    var='var',
    uns='uns'
)

adata_normalize_matrix = ecdf_qnorm_genes(adata, batch_key=par['batch_key'])

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
    X=adata_normalize_matrix,
    obs=adata.obs.copy(),
    var=adata.var.copy(),
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
    }
)

output.write_h5ad(par['output'], compression='gzip')
