import os
import re
import pandas as pd
from datetime import date
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List


def get_date(neuron: str) -> date:
    """
    Extracts date data from singular index
    :param neuron: string representing single neuron
    :return: date data
    """
    date_pattern = re.compile(r'(20\d\d)(\d\d)(\d\d)')
    match = date_pattern.search(neuron)
    if not match:
        raise ValueError(f'Date not found in {neuron}')
    year, month, day = tuple([int(x) for x in match.groups()])
    if month == 0:
        month = 1
    return date(year, month, day)


def get_exp(neuron: str) -> str:
    """
    Get experiment (DD || LD) applied to the neuron
    :param neuron:
    :return: string 'DD' or 'LD'
    """
    exp_re = re.compile(r"(_)(DD|LD)(_)")
    match = exp_re.search(neuron)
    if not match:
        raise ValueError(f'Experiment not found in {neuron}')
    return match.group(2)


def get_time(neuron: str) -> str:
    """
    Same as above, capturing ZT or CT time
    :param neuron: single neuron name
    :return: ZTXX or CTXX
    """
    time_pattern = re.compile(r'(CT|ZT)(\d\d)')
    match = time_pattern.search(neuron)
    if not match:
        raise ValueError(f'Time not found in {neuron}')
    return match.group()


def read_genexp_files(genes: List[str],
                      DATA_PATH: str = r'../dataset/',
                      ANNOT_PATH: str = 'clock_neuron_clusters.csv') -> pd.DataFrame:
    """
    Reads through all .csv files in DATA_PATH 

    Args:
        genes (List[str]): relevant genes to look for
        DATA_PATH (str, optional): directory with available .csv files. Defaults to r'../dataset/'.
        ANNOT_PATH (str, optional): .csv file with further annotations. Defaults to 'clock_neuron_clusters.csv'.

    Returns:
        pd.DataFrame: each row is a single cell, columns indicate either gene expression or annotations
    """

    annot_df = pd.read_csv(ANNOT_PATH, index_col=0)
    mapper = lambda x : x[1:] # removes leading 'x' char in idx strings
    genexp_df = None

    # Part 1: iterate over all files in dir, and get relevant data
    print('Going through files...')
    for filename in os.listdir(DATA_PATH):
        new_df = pd.read_csv(os.path.join(DATA_PATH, filename), index_col=0)
        new_df = new_df.loc[genes]\
                        .rename(mapper, axis='columns')\
                        .T\
                        .merge(annot_df[['Idents']], left_index=True, right_index=True, how='left')\
                        .dropna(subset=['Idents'])
            
        if genexp_df is None:
            genexp_df = new_df
            continue
        genexp_df = pd.concat([genexp_df, new_df])
    
    # Part 2: add further annotations
    print('Getting annotations...')
    indeces = genexp_df.index.to_series()
    genexp_df['date'] = indeces.map(get_date)
    genexp_df['experiment'] = indeces.map(get_exp)
    genexp_df['exp_time'] = indeces.map(get_time)

    return genexp_df


def df_to_anndata(df: pd.DataFrame) -> AnnData:
    """
    Transform single cell genetic expression dataframe to
    the annotated data class, by separating gene columns from 
    annotation fields

    Args:
        df (pd.DataFrame): indeces are single cells, columns are either genes or annot

    Returns:
        AnnData: same info. as the df, but now annotations are stored separetely
    """
    df_cols = df.columns
    annot_fields = ['Idents', 'date', 'experiment', 'exp_time']
    assert any([field in df_cols for field in annot_fields])
    annot_cols = [col for col in annot_fields if col in df.columns]
    
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
    # genes = ["DIP-gamma", "DIP-beta", "DIP-delta", "DIP-theta"]
    # genexp_df = read_genexp_files(genes)
    # print(genexp_df.shape)
    # print(genexp_df.head(10))
    
    # OUT_DIR = 'data_subsets/'
    # save = input("Save to csv? [y/n]: ")
    # if save.strip().lower() == 'y':
    #     OUT_NAME = input('File name: ')
    #     if '.csv' not in OUT_NAME:
    #         OUT_NAME += '.csv'
    #     genexp_df.to_csv(os.path.join(OUT_DIR, OUT_NAME))

    # Test 2: converto to anndata
    read_df = pd.read_csv(r"data_subsets\DIP-genes.csv", index_col=0)
    genex_ad = df_to_anndata(read_df)
    print(genex_ad)


if __name__ == '__main__':
    main()