import os
import re
import pandas as pd
from datetime import date
from anndata import AnnData
import time
from typing import List


with open('all_genes.txt', 'r') as f:
    ALL_GENES = f.read().splitlines()[1:]


def read_genexp_files(genes: List[str] | None = None,
                      data_path: str = r'../dataset/') -> pd.DataFrame:
    """
    Reads through all .csv files in DATA_PATH 

    Args:
        genes (List[str]): relevant genes to look for
        data_path (str, optional): directory with available .csv files. Defaults to r'../dataset/'.

    Returns:
        pd.DataFrame: each row is a single cell, columns indicate either gene expression or annotations
    """
    if genes is None:
        genes = ['amon']

    # Part 1: iterate over all files in dir, and get relevant data
    genexp_df = None
    print('Going through files...')
    for filename in os.listdir(data_path):
        new_df = pd.read_csv(os.path.join(data_path, filename), index_col=0, keep_default_na=False)
        new_df = new_df.rename(lambda x: x[1:], axis='columns').T  # removes leading 'x' char in idx strings
        try:
            new_df = new_df[genes]
        except KeyError:
            continue
        if genexp_df is None:
            genexp_df = new_df
            continue
        genexp_df = pd.concat([genexp_df, new_df])

    # Part 2: add further annotations
    return genexp_df


def df_to_anndata(df: pd.DataFrame, annot_cols: List[str] | None = None) -> AnnData:
    """
    Transform single cell genetic expression dataframe to
    the annotated data class, by separating gene columns from 
    annotation fields

    Args:
        df (pd.DataFrame): indeces are single cells, columns are either genes or annot
        annot_cols: list of annotation columns

    Returns:
        AnnData: same info. as the df, but now annotations are stored separetely
    """
    if annot_cols is None:
        annot_cols = ['experiment', 'Repeats', 'condition', 'date', 'time', 'Idents']

    assert all([field in df.columns for field in annot_cols])

    obs = df[annot_cols]
    df = df.drop(annot_cols, axis='columns')
    adata = AnnData(X=df, obs=obs, var=df.columns.to_frame())
    return adata


def anndata_to_df(adata: AnnData) -> pd.DataFrame:
    """
    Convert back to df from anndata without losing the observation columns
    Args:
        adata (AnnData): annotated data object, with 'vars' and 'obs'

    Returns:
        pd.DataFrame: equivalent dataframe
    """
    expr_df = adata.to_df()
    obs_df = adata.obs
    combined_df = pd.concat([expr_df, obs_df], axis=1)
    return combined_df


def main():
    # Test 1: read from raw data
    genes = ["DIP-gamma", "DIP-beta", "DIP-delta", "DIP-theta", "dpr8"]
    annot_path = 'neuron_annotations.csv'
    annot_df = pd.read_csv(annot_path, index_col=0)
    start = time.time()
    genexp_df = read_genexp_files(genes)
    end = time.time()
    print(f'Read genexp took {end - start:.2f} seconds')
    full_df = genexp_df.join(annot_df, how='right')
    print(full_df.shape)
    print(full_df.head(10))

    OUT_DIR = 'data_subsets/'
    save = input("Save to csv? [y/n]: ")
    if save.strip().lower() == 'y':
        OUT_NAME = input('File name: ')
        if '.csv' not in OUT_NAME:
            OUT_NAME += '.csv'
        genexp_df.to_csv(os.path.join(OUT_DIR, OUT_NAME))

    # Test 2: convert to anndata
    genexp_ad = df_to_anndata(full_df)
    print(genexp_ad)


if __name__ == '__main__':
    main()
