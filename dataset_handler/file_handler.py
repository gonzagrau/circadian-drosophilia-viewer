import os
import gdown
import pandas as pd
from anndata import AnnData
import time
from typing import List, Tuple
from scanpy import read_h5ad


DD_DATA_URL = r'https://drive.google.com/uc?id=1jewJcYuaPyVE327VI0viPP-9kjtFCx_7'
LD_DATA_URL = r' https://drive.google.com/uc?id=1KF4ce-rdJ7-PPumVlka2bFwoxV0Db1av'


def load_h5ad_files(datadir: str = 'dataset',
                    DD_url: str = DD_DATA_URL,
                    LD_url: str = LD_DATA_URL) -> Tuple[AnnData, AnnData]:
    """
    Fetches data from h5ad files. If they do not exist, they are downloaded.
    Args:
        datadir: where to find the h5ad files, if they exist
        DD_url: url to download the data from dark-dark experiments
        LD_url: url to download the data from light-dark experiments

    Returns:
        a tuple containing anndata objects of LD and DD experiments, respectively.
    """
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    datapaths = [os.path.join(datadir, f'dataset_{cond}.h5ad') for cond in ['DD', 'LD']]
    if not all([os.path.exists(dpath) for dpath in datapaths]):
        gdown.download(DD_url, datapaths[0], quiet=False)
        gdown.download(LD_url, datapaths[1], quiet=False)

    ad_LD = read_h5ad(os.path.join(datadir, 'dataset_LD.h5ad'))
    ad_DD = read_h5ad(os.path.join(datadir, 'dataset_DD.h5ad'))

    return ad_LD, ad_DD

def read_genexp_files(genes: List[str] | None = None,
                      data_path: str = r'dataset/') -> pd.DataFrame:
    """
    Reads through all .csv files in DATA_PATH 

    Args:
        genes (List[str]): relevant genes to look for
        data_path (str, optional): directory with available .csv files. Defaults to r'dataset/'.

    Returns:
        pd.DataFrame: each row is a single cell, columns indicate either gene expression or annotations
    """
    if genes is None:
        genes = ['amon']

    # Part 1: iterate over all files in dir, and get relevant data
    genexp_df = None
    print('Going through files...')
    for filename in os.listdir(data_path):
        new_df = pd.read_csv(os.path.join(data_path, filename), index_col=0) #, keep_default_na=False)
        new_df = new_df.rename(lambda x: x[1:], axis='columns').T  # removes leading 'x' char in idx strings
        try:
            new_df = new_df[genes]
        except KeyError:
            continue
        if genexp_df is None:
            genexp_df = new_df
            continue
        genexp_df = pd.concat([genexp_df, new_df])

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
    annot_path = '../neuron_annotations.csv'
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
    #main()
    load_h5ad_files()
