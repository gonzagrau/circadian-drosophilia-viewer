import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import scanpy as sc
from scanpy import AnnData
import anndata as ad
from dataset_handler import anndata_to_df


sns.set_palette('deep')
# mpl.use('TkAgg')

def preprocess_pipeline(adata: AnnData,
                        patterns_to_exclude: list | None = None,
                        exclude_highly_expressed: bool = False,
                        log_normalize: bool = False) -> AnnData:
    """
    Excludes genes based on certain name patterns
    Args:
        adata: gene expression annotated data object
        patterns_to_exclude: list of patterns to exclude. if None, defaults to a pre-specified pattern defined below
        exclude_highly_expressed: kwarg for normalize_total. defaults to False
        log_normalize: second normalization, optional, applying a log transform. defaults to False.
    Returns:
        modified, normalized anndata object. no changes are made inplace to the original.
    """
    if patterns_to_exclude is None:
        patterns_to_exclude = ['mt', 'rpls', 'rRNA', 'tRNA', 'ERCC', 'EGFP']

    # Filter genes
    gene_names = adata.var_names.astype(str)
    filtered_genes = ~gene_names.str.contains('|'.join(patterns_to_exclude), regex=True)
    adata = adata[:, filtered_genes].copy()

    # Normalize per cell (same total count for each cell)
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=exclude_highly_expressed)

    # Log normalize
    if log_normalize:
        sc.pp.log1p(adata)

    return adata


def test_preprocess_pipeline(adata: AnnData) -> None:
    """
    Makes some basic plots using a preprocessed anndata object.
    Args:
        adata: preprocessed anndata object.

    Returns:
        None
    """

    # Seleccionar los genes de interes
    genes_of_interest = ['DIP-gamma', 'DIP-beta', 'DIP-delta', 'DIP-theta', 'dpr8']

    # Crear el dotplot
    sc.pl.dotplot(adata, var_names=genes_of_interest, groupby='Idents', standard_scale='var',)

    # Filtrar las celulas que pertenecen al cluster 'sLNVs'
    adata_slNvs = adata[adata.obs['Idents'] == '2:s_LNv', :]
    # adata_llNvs = adata[adata.obs['Idents'] == '25:l_LNv', :]

    # Obtener los datos de expresion del gen
    df = anndata_to_df(adata_slNvs[:, adata_slNvs.var['gene_names'].isin(genes_of_interest)])

    # Crear el grafico usando Seaborn pointplot
    plt.figure(figsize=(8, 5))
    sns.pointplot(data=df, x='time', y='DIP-beta', estimator='mean', errorbar='se')
    plt.xlabel('Time (ZT)')
    plt.ylabel('Mean Expression of DIP-beta')
    plt.title('Mean Expression in sLNVs Cluster Over Time')
    plt.grid(True)
    plt.show()


def main():
    adata_ld = sc.read_h5ad('dataset_LD.h5ad')
    adata_dd = sc.read_h5ad('dataset_DD.h5ad')
    adata = ad.concat([adata_ld, adata_dd], join='inner')
    genes_of_interest = ['DIP-gamma', 'DIP-beta', 'DIP-delta', 'DIP-theta', 'dpr8']
    adata = adata[:, adata.var.index.isin(genes_of_interest)]
    return
    filtered_adata = preprocess_pipeline(adata)
    test_preprocess_pipeline(filtered_adata)


if __name__ == '__main__':
    main()
