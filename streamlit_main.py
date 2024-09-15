import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import List, Tuple
from copy import deepcopy
from dataset_handler import df_to_anndata, anndata_to_df, load_h5ad_files
from pipeline import preprocess_pipeline
from random import randint
from datetime import datetime

# mpl.use('TKAgg')
sns.set_palette('deep')
heatmap_palette = sns.diverging_palette(240, 50, l=30, as_cmap=True)
ANNOT_PATH = r"neuron_annotations.csv"

@st.cache_data
def get_full_adata() -> ad.AnnData:
    """
    Get the anndata objects from the h5ad files in the dataset folder, and caches it.
    Refer to dataset_handler.load_h5ad_files() for further information
    """
    ad_LD, ad_DD = load_h5ad_files()
    adata = ad.concat([ad_LD, ad_DD], join='inner')
    filtered_adata = preprocess_pipeline(adata)
    st.session_state['full_adata'] = filtered_adata
    return filtered_adata

@st.cache_data
def fetch_data(choice: List[str]) -> None:
    """
    Reads data from csv files and stores it in session state
    Args:
        choice (List[str]): choice of genes to select from dataset
    """
    print(f"User gene choice: {choice}")
    st.session_state['genes'] = choice
    with st.spinner(text='Fetching data...'):
        filtered_adata = get_full_adata()
        final_adata = filtered_adata[:, filtered_adata.var_names.isin(choice)].copy()
        # cache variables
        st.session_state['adata'] = final_adata
        st.session_state['dataframe'] = anndata_to_df(final_adata)

    # Save all found clusters
    idents = [i for i in st.session_state.dataframe['Idents'].unique()]
    idents.sort(key=lambda x: int(x.split(':')[0]))
    st.session_state['Idents'] = idents


def idents_multiselect(key: str) -> List[str]:
    """
    Widget to select one or more neuron clusters
    Returns:
        List[str]: selected clusters
    """
    idents = st.session_state['Idents']
    return st.multiselect("Select a cluster", idents, default=idents[0], key=key)


def make_dotplots() -> None:
    """
    Interactive dotplot builder for gene expression
    """
    rightnow = datetime.now().microsecond
    keynum = rightnow*randint(0, rightnow)

    # Backend data reordering
    adata = deepcopy(st.session_state['adata'])
    exps = ['LD', 'DD']

    st.write('## Select data to display')
    condition_choice = st.multiselect("Select condition", exps, default='LD', key=f"cond_{keynum}")
    adata_plot = adata[adata.obs['condition'].isin(condition_choice)].copy()

    # Group by
    group = st.radio('Group by:', ['time', 'Idents'], key=f'groupby_{keynum}')
    extra_str = ''
    if group == 'time':
        idents = st.multiselect("Select a cluster",
                                st.session_state['Idents'],
                                default=st.session_state['Idents'][0],
                                key=f"cluster_{keynum}")
        adata_plot = adata_plot[adata_plot.obs['Idents'].isin(idents)].copy()
        extra_str = f" for {','.join(idents)}"

    # Choose gene subset
    var_names = st.multiselect("Select genes in order", adata_plot.var_names.unique(),
                               default=st.session_state['genes'],
                               key=f"vnames_{keynum}")

    # Plot!
    swap_axes = st.checkbox('Swap axes')
    title = f"Gene expression by {group} at {','.join(condition_choice)}" + extra_str
    dotplot = sc.pl.DotPlot(adata_plot,
                            var_names=var_names,
                            groupby=group,
                            standard_scale='var',
                            vmin=-1,
                            vmax=2,
                            var_group_rotation=0.,
                            edgecolors=None,
                            mean_only_expressed=True,
                            title=title,
                            cmap='Reds',
                            linewidth=0.)
    if swap_axes:
        dotplot.swap_axes()
    dotplot.make_figure()
    st.session_state['dotplot'] = dotplot.fig


def make_pointplots() -> None:
    """
    Makes pointplots of expression vs. time for every gene
    """
    adata = deepcopy(st.session_state['adata'])

    # First, pick clusters
    st.write('## Design your hourly expression pointplots')
    id_choice = st.multiselect("Pick a cluster", 
                                st.session_state['Idents'], 
                                default=st.session_state['Idents'][0])

    df = anndata_to_df(adata[adata.obs['Idents'].isin(id_choice)].copy())
    df['time'] = df['time'].apply(lambda x: x[2:])  # (ZT||CT)XX -> XX

    # Plot
    figures = []
    gene_palette = {"LD": 'turquoise', "DD": 'gray'}
    for i, gene in enumerate(adata.var_names):
        fig, ax = plt.subplots()
        sns.pointplot(df,
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


def make_heterogeneity_heatmap() -> None:
    """
    Designs a heatmap and saves it to the session state
    """
    st.write('Inner cluster heterogeneity at a given time')
    adata = deepcopy(st.session_state['adata'])
    take_log = st.toggle('Apply logarithm')
    df = anndata_to_df(adata)
    
    # Pick clusters and times
    times = [t for t in df['time'].unique()]
    id_choice = st.multiselect("Choose a cluster", 
                               st.session_state['Idents'], 
                               default=st.session_state['Idents'][1])
    t_choice = st.selectbox('Pick the time', times, placeholder='ZT or CT...')
    
    # Select relevant data
    df = df[(df["Idents"].isin(id_choice)) & (df['time'] == t_choice)]
    df = df.select_dtypes(include=('float', 'int'))
    df = df - df.mean(axis=0)

    # Plot
    fig, ax = plt.subplots()
    ax.set_title(f"Heterogeneity for {', '.join(id_choice)} at {t_choice}")
    heatmap = sns.heatmap(df,
                          ax=ax,
                          cmap=heatmap_palette,
                          cbar=True,
                          center=0,
                          vmin=-3,
                          vmax=3,
                          yticklabels=False, 
                          xticklabels=True)
    st.session_state['heatmap'] = heatmap.figure   


def make_matrixplots() -> None:
    """
    Design matrix plots (i.e: heatmaps) just like dotplots
    """

    # Backend data reordering
    adata = deepcopy(st.session_state['adata'])
    exps = ['LD', 'DD']

    st.write('## Select data to display')
    condition_choice = st.multiselect("Select condition", exps, default='DD')
    adata_plot = adata[adata.obs['condition'].isin(condition_choice)].copy()

    # Group by
    group = st.radio('Group by:', ['time', 'Idents'])
    extra_str = ''
    if group == 'time':
        idents = st.multiselect("Select one or more clusters",
                                st.session_state['Idents'],
                                default=st.session_state['Idents'][0])
        adata_plot = adata_plot[adata_plot.obs['Idents'].isin(idents)].copy()
        extra_str = f" for {','.join(idents)}"

    # Choose gene subset
    var_names = st.multiselect("Determine gene order", adata_plot.var_names.unique(),
                               default=st.session_state['genes'])

    swap_axes = st.checkbox('Show groups in X axis')

    # Plot!
    title = f"Mean gene expression by {group} at {','.join(condition_choice)}" + extra_str
    matrixplot = sc.pl.MatrixPlot(adata_plot,
                                    var_names=var_names,
                                    groupby=group,
                                    var_group_rotation=0.,
                                    title=title,
                                    vmin=-1,
                                    vmax=2,
                                    cmap='viridis')
    if swap_axes:
        matrixplot.swap_axes()
    matrixplot.make_figure()
    st.session_state['matrixplot'] = matrixplot.fig


def main():
    # Page Title
    st.set_page_config(page_title="CircDrosView", page_icon="bar-chart")
    st.write("""
    # Interactive visualizer for *A transcriptomic taxonomy of* Drosophila *circadian neurons around the clock*
             """)

    # Initialization
    with st.spinner(text='Initializing variables...'):
        with open('all_genes.txt', 'r') as f:
            genes = f.read().splitlines()
        if 'genes' not in st.session_state:
            st.session_state['genes'] = []
            st.session_state['full_adata'] = sc.AnnData()
            st.session_state['dataframe'] = pd.DataFrame()
            st.session_state['adata'] = sc.AnnData()
            st.session_state['Idents'] = []
            st.session_state['dotplot'] = None
            st.session_state['pointplots'] = []
            st.session_state['heatmap'] = None
            st.session_state['matrixplot'] = None
        get_full_adata()

    # Gene selection
    st.write('## Step 1: select genes to analyze')
    with st.form(key="input_parameters"):
        choice = st.multiselect("Choose your genes", genes)
        submit = st.form_submit_button("Submit")

    if submit:
        fetch_data(choice)

    if not st.session_state['dataframe'].empty:
        st.write("Fetched Data:")
        df = st.session_state['dataframe']
        st.dataframe(df)

        tab_dot, tab_point, tab_het, tab_mat = st.tabs(['Dot plots', 'Point plots', 'Heterogeneity heatmap', 'Matrix plot'])
        # Dotplots
        with tab_dot:
            plt.close()
            make_dotplots()
            if st.session_state['dotplot'] is not None:
                st.pyplot(st.session_state['dotplot'])
            
        # Pointplots
        with tab_point:
            plt.close()
            make_pointplots()
            if len(st.session_state['pointplots']):
                for figure in st.session_state['pointplots']:
                    plt.legend()
                    st.pyplot(figure)

        # Heterogeneity Heatmaps
        with tab_het:
            plt.close()
            make_heterogeneity_heatmap()
            if st.session_state['heatmap'] is not None:
                st.pyplot(st.session_state['heatmap'])

        with tab_mat:
            plt.close()
            make_matrixplots()
            if st.session_state['matrixplot'] is not None:
                st.pyplot(st.session_state['matrixplot'])
    else:
        st.write("Please select data to fetch.")


if __name__ == '__main__':
    main()
