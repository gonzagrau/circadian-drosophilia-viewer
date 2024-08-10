import streamlit as st
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import scanpy as sc
from dataset_handler import read_genexp_files, df_to_anndata, anndata_to_df
from typing import List

mpl.use('agg')
sns.set_palette('deep')
ANNOT_PATH = r"neuron_annotations.csv"


def fetch_data(choice: List[str]) -> None:
    """
    Reads data from csv files and stores it in session state
    Args:
        choice (List[str]): choice of genes to select from dataset
    """
    print(choice)
    with st.spinner(text='Fetching data...'):
        annot_df = pd.read_csv(ANNOT_PATH, index_col=0)
        data_df = read_genexp_files(choice)
        df = data_df.join(annot_df, how='right')
        st.success('Done')
    print(df.shape)
    st.session_state['data'] = df


def make_dotplots(df: pd.DataFrame) -> None:
    """
    Interactive dotplot builder
    """
    adata = df_to_anndata(df)
    idents = [ident for ident in adata.obs['Idents'].unique()]
    exps = ['LD', 'DD']
    idents.sort(key=lambda x: int(x.split(':')[0]))

    st.write('## Step 2: design your dotplots')
    condition_choice = st.multiselect("Select condition condition", exps, default='LD')
    df_exp = df[df['condition'].isin(condition_choice)]
    group = st.radio('Group by:', ['time', 'Idents'])

    if group == 'Idents':
        adata_dotplot = df_to_anndata(df_exp)
        extra_str = ''
        swap = False

    else:  # elif group == 'time':
        idents = st.multiselect("Select a cluster", idents, default=idents[0])
        df_group = df_exp[(df_exp['Idents'].isin(idents))]
        adata_dotplot = df_to_anndata(df_group)
        extra_str = f"for {','.join(idents)}"
        swap = True

    title = f"Gene expression by {group} at {','.join(condition_choice)}" + extra_str
    dotplot = sc.pl.DotPlot(adata_dotplot,
                            var_names=adata.var_names,
                            groupby=group,
                            standard_scale='var',
                            vmin=-1,
                            vmax=1,
                            var_group_rotation=0.,
                            edgecolors=None,
                            mean_only_expressed=True,
                            title=title,
                            cmap='Reds',
                            linewidth=0.)
    if swap:
        dotplot.swap_axes()
    dotplot.make_figure()
    st.session_state['dotplot'] = dotplot.fig


def make_pointplots(adata: sc.AnnData) -> None:
    """
    Makes pointplots for every gene
    :param adata: contains gene expression and annotations
    :return:
    """
    idents = [ident for ident in adata.obs['Idents'].unique()]
    idents.sort(key=lambda x: int(x.split(':')[0]))

    # First, pick clusters
    st.write('## Step 2: design your pointplots')
    id_choice = st.multiselect("Select clusters", idents, default=idents[0])

    # Then, preprocess data
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
    df = anndata_to_df(adata)
    df['time'] = df['time'].apply(lambda x: x[2:])  # (ZT||CT)XX -> XX
    df_clust = df[df['Idents'].isin(id_choice)]

    # Plot
    figures = []
    gene_palette = {"LD": 'turquoise', "DD": 'gray'}
    for i, gene in enumerate(adata.var_names):
        fig, ax = plt.subplots()
        sns.pointplot(df_clust,
                      x='time',
                      y=gene,
                      hue='condition',
                      estimator='mean',
                      errorbar='se',
                      palette=gene_palette,
                      capsize=0.2,
                      linewidth=1.5,
                      ax=ax,)
        ax.set_ylabel('Gene expression (TP10K)')
        ax.set_title(f"{gene} expression in {', '.join(id_choice)}")
        figures.append(fig)

    # Save
    st.session_state['pointplots'] = figures


def main():
    # Page Title
    st.set_page_config(page_title="CircDrosView", page_icon="bar-chart")
    st.write("""
    # Interactive visualizer for *A transcriptomic taxonomy of* Drosophilia *circadian neurons around the clock*
             """)

    # Initialization
    with st.spinner(text='Initializing variables...'):
        with open('all_genes.txt', 'r') as f:
            genes = f.read().splitlines()
        if 'data' not in st.session_state:
            st.session_state['data'] = pd.DataFrame()
            st.session_state['dotplot'] = None
            st.session_state['pointplots'] = []

    # Gene selection
    st.write('## Step 1: select genes to analyze')
    with st.form(key="input_parameters"):
        choice = st.multiselect("Choose your genes", genes)
        submit = st.form_submit_button("Submit")

    if submit:
        fetch_data(choice)

    if not st.session_state['data'].empty:
        st.write("Fetched Data:")
        df = st.session_state['data']
        st.dataframe(df)

        # Dotplots
        make_dotplots(df)
        if st.session_state['dotplot'] is not None:
            st.pyplot(st.session_state['dotplot'])
        
        # Pointplots
        adata = df_to_anndata(df)
        make_pointplots(adata)

        if len(st.session_state['pointplots']):
            for figure in st.session_state['pointplots']:
                st.pyplot(figure)

    else:
        st.write("Please select data to fetch.")


if __name__ == '__main__':
    main()
