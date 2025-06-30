from typing import Any, Callable, Optional, Dict

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


def ecdf_qnorm_samples(
    adata: sc.AnnData,
    layer: Optional[str] = None,
    summary_function: Callable = np.nanmean,
    merge_grid_size: int = 1000,
    verbose: bool = True,
):
    """
    Run quantile normalization on samples
    :param adata: Input adata
    :param layer: layer name
    :param summary_function: function to norm with, such as np.nanmean np.nanmedian, etc
    :param merge_grid_size:
    :param verbose: Whether to print progress
    :return: a csr_matrix with integrated values
    """
    if layer:
        data_in = adata.layers[layer]
    else:
        data_in = adata.X

    num_samples, num_genes = data_in.shape
    exp_output = sparse.lil_matrix(np.zeros((num_samples, num_genes)))
    merge_grid = np.linspace(0, 1, merge_grid_size + 1)
    freq_matrix = sparse.lil_matrix(np.zeros((num_samples, num_genes)))

    merge_grid_matrix = np.zeros((num_samples, len(merge_grid)))
    data_in = data_in.tocsr() if sparse.issparse(data_in) else data_in

    with tqdm(total=num_samples, disable=not verbose) as pbar:
        for sample_index in range(num_samples):
            sample_data = (
                data_in.getrow(sample_index).toarray().flatten()
                if sparse.issparse(data_in)
                else data_in[sample_index, :]
            )

            ranked_data = rankdata(sample_data, method="average")
            normalized_ranks = ranked_data - np.min(ranked_data)
            normalized_ranks = normalized_ranks / np.maximum(
                np.max(normalized_ranks), 1
            )

            unique_sample_data, unique_indices = np.unique(
                sample_data, return_index=True
            )
            x_normalized_ranks = normalized_ranks[unique_indices]
            x_sample_data = sample_data[unique_indices]

            merge_grid_matrix[sample_index, :] = np.maximum(
                np.interp(
                    merge_grid,
                    x_normalized_ranks,
                    x_sample_data,
                    left=np.nan,
                    right=np.nan,
                ),
                0,
            )

            freq_matrix[sample_index, :] = normalized_ranks
            pbar.update(1)

    merge_grid_matrix_summary = summary_function(merge_grid_matrix, axis=0)
    # Generate normalized matrix
    non_zero_indices = freq_matrix.nonzero()
    exp_output[non_zero_indices] = np.maximum(
        np.interp(
            freq_matrix[non_zero_indices].toarray(),
            merge_grid,
            merge_grid_matrix_summary,
            left=np.nan,
            right=np.nan,
        ),
        0,
    )
    return exp_output.tocsr()
