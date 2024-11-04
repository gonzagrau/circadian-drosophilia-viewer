import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import List
from copy import deepcopy
from dataset_handler.file_handler import anndata_to_df, load_filtered_adata

####################################################################################
#                            GLOBAL VARIABLES                                      #
####################################################################################
global COND_PALETTE


def color_theme_setup() -> None:
    """
    Change default plt config to match selected colorthem
    """

    global COND_PALETTE
    # Detect the user's theme
    is_dark_mode = st.get_option("theme.base") == "dark"

    # Set the color palette and matplotlib rcParams based on the theme
    if is_dark_mode:
        sns.set_palette('dark')
        COND_PALETTE = {"LD": 'lightblue', "DD": 'gray'}
        plt.rcParams.update({
            "axes.facecolor": "black",
            "axes.edgecolor": "white",
            "axes.labelcolor": "white",
            "figure.facecolor": "black",
            "figure.edgecolor": "black",
            "grid.color": "gray",
            "text.color": "white",
            "xtick.color": "white",
            "ytick.color": "white",
            "savefig.facecolor": "black",
            "savefig.edgecolor": "black"
        })
    else:
        sns.set_palette('deep')
        COND_PALETTE = {"LD": 'turquoise', "DD": 'gray'}
        plt.rcParams.update({
            "axes.facecolor": "white",
            "axes.edgecolor": "black",
            "axes.labelcolor": "black",
            "figure.facecolor": "white",
            "figure.edgecolor": "white",
            "grid.color": "gray",
            "text.color": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white"
        })


# List of colormaps
COLORMAPS = plt.colormaps()


####################################################################################
#                            CACHED FUNCTIONS                                      #
####################################################################################

@st.cache_data
def get_full_adata() -> ad.AnnData:
    """
    Get the anndata objects from the h5ad files in the dataset folder, and caches it.
    Refer to dataset_handler.load_h5ad_files() for further information
    """
    adata = load_filtered_adata()
    st.session_state['full_adata'] = adata
    return adata


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


@st.cache_data
def compile_scplot(plot_type: str,
                   _adata_plot: ad.AnnData,
                   group: str,
                   condition_choice: List[str],
                   var_names: List[str],
                   extra_str: str,
                   swap_axes: bool,
                   cmap: str = 'Reds') -> plt.figure:
    """
    Given a predesigned plot object, make its figure and save it
    """
    assert plot_type in ['dotplot', 'matrixplot']
    title = f"Gene expression by {group} at {', '.join(condition_choice)}" + extra_str

    if plot_type == 'dotplot':
        plot_obj = sc.pl.DotPlot(_adata_plot,
                                 var_names=var_names,
                                 groupby=group,
                                 standard_scale='var',
                                 vmin=-1,
                                 vmax=2,
                                 var_group_rotation=0.,
                                 edgecolors=None,
                                 mean_only_expressed=True,
                                 title=title,
                                 cmap=cmap,
                                 linewidth=0.)
    elif plot_type == 'matrixplot':
        plot_obj = sc.pl.MatrixPlot(_adata_plot,
                                    var_names=var_names,
                                    groupby=group,
                                    standard_scale='var',
                                    var_group_rotation=0.,
                                    title=title,
                                    vmin=-1,
                                    vmax=2,
                                    cmap=cmap)
    else:
        raise ValueError("Invalid plot type")

    if swap_axes:
        plot_obj.swap_axes()
    plot_obj.make_figure()
    return plot_obj.fig


@st.cache_data
def compile_pointplots(df: pd.DataFrame,
                       genes: List[str],
                       id_choice: List[str]) -> List[plt.figure]:
    """
    Generates day-long expression plots for all selected genes in both conditions.
    """
    figures = []
    for i, gene in enumerate(genes):
        fig, ax = plt.subplots()
        sns.pointplot(df,
                      x='time',
                      y=gene,
                      hue='condition',
                      estimator='mean',
                      errorbar='se',
                      palette=COND_PALETTE,
                      capsize=0.2,
                      linewidth=1.5,
                      ax=ax, )
        ax.set_ylabel('Gene expression (TP10K)')
        ax.set_title(f"{gene} expression in {', '.join(id_choice)}")
        figures.append(fig)
    return figures


@st.cache_data
def compile_heterogeneity(df: pd.DataFrame,
                          id_choice: List[str],
                          t_choice: List[str],
                          cmap: str = 'viridis') -> plt.figure:
    """
    Generates a heatmap of heterogeneity for a given cluster at a given time
    Args:
        df: preprocessed dataframe for gene expression
        id_choice: chosen cluster
        t_choice: chosen time
        cmap: colormap to use

    Returns:
        plt.figure: the heatmap figure
    """

    fig, ax = plt.subplots()
    ax.set_title(f"Heterogeneity for {', '.join(id_choice)} at {t_choice}")
    hmap = sns.heatmap(df,
                       ax=ax,
                       cmap=cmap,
                       cbar=True,
                       center=0,
                       vmin=-3,
                       vmax=3,
                       yticklabels=False,
                       xticklabels=True)
    return hmap.figure


@st.cache_data
def compile_tsne_plot(cmap: str = 'viridis') -> None:
    """
    Makes t-SNE plots for the selected data
    Returns: None, but stores tsne plots as figures to the cache
    """
    # General t-SNE plot for all idents
    fulladata = deepcopy(st.session_state['full_adata'])
    ident_tsne = sc.pl.tsne(fulladata, color='Idents', return_fig=True)
    st.session_state['tsne_plots'] = [ident_tsne]

    for gene in st.session_state['genes']:
        gene_tsne = sc.pl.tsne(fulladata, color=gene, cmap=cmap, return_fig=True)
        st.session_state['tsne_plots'].append(gene_tsne)


@st.cache_data
def compile_umap_plot(cmap: str = 'viridis') -> None:
    """
    Makes UMAP plots for the selected data
    Returns: None, but stores umap plots as figures to the cache
    """
    # General UMAP plot for all idents
    fulladata = deepcopy(st.session_state['full_adata'])
    ident_umap = sc.pl.umap(fulladata, color='Idents', return_fig=True)
    st.session_state['umap_plots'] = [ident_umap]

    for gene in st.session_state['genes']:
        gene_umap = sc.pl.umap(fulladata, color=gene, cmap=cmap, return_fig=True)
        st.session_state['umap_plots'].append(gene_umap)


####################################################################################
#                        FUNCTIONS INVOLVING WIDGETS                               #
####################################################################################

def make_dotplots() -> None:
    """
    Interactive dotplot builder for gene expression
    """
    # Backend data reordering
    adata = deepcopy(st.session_state['adata'])
    exps = ['LD', 'DD']

    st.write('## Select data to display')
    condition_choice = st.multiselect("Select light condition", exps, default='LD', key="cond_dot")
    adata_plot = adata[adata.obs['condition'].isin(condition_choice)].copy()

    # Group by
    group = st.radio('Group by:', ['time', 'Idents'], key='groupby_dot')
    extra_str = ''
    if group == 'time':
        idents = st.multiselect("Select a cluster",
                                st.session_state['Idents'],
                                default=st.session_state['Idents'][0],
                                key="cluster_dot")
        adata_plot = adata_plot[adata_plot.obs['Idents'].isin(idents)].copy()
        extra_str = f" for {','.join(idents)}"

    # Choose gene subset
    var_names = st.multiselect("Select genes in order", adata_plot.var_names.unique(),
                               default=st.session_state['genes'],
                               key="vnames_dot")

    # Select colormap
    cmap = st.selectbox("Select colormap", COLORMAPS, key='cmap_dot')

    # Plot!
    swap_axes = st.checkbox('Swap axes', key='swap_dot')
    fig = compile_scplot('dotplot',
                         _adata_plot=adata_plot,
                         group=group,
                         condition_choice=condition_choice,
                         var_names=var_names,
                         extra_str=extra_str,
                         swap_axes=swap_axes,
                         cmap=cmap)
    st.session_state['dotplot'] = fig


def make_pointplots() -> None:
    """
    Makes pointplots of expression vs. time for every gene
    """
    adata = deepcopy(st.session_state['adata'])

    # First, pick clusters
    st.write('## Design your hourly expression pointplots')
    id_choice = st.multiselect("Pick a cluster",
                               st.session_state['Idents'],
                               default=st.session_state['Idents'][0],
                               key="cluster_point")

    df = anndata_to_df(adata[adata.obs['Idents'].isin(id_choice)].copy())
    df['time'] = df['time'].apply(lambda x: x[2:])  # (ZT||CT)XX -> XX

    # Plot
    figures = compile_pointplots(df, list(adata.var_names), id_choice)

    # Save
    st.session_state['pointplots'] = figures


def make_heterogeneity_heatmap() -> None:
    """
    Designs a heatmap and saves it to the session state
    """

    st.write('Inner cluster heterogeneity at a given time')
    adata = deepcopy(st.session_state['adata'])
    take_log = st.toggle('Apply logarithm')
    if take_log:
        sc.pp.log1p(adata)
    df = anndata_to_df(adata)

    # Pick clusters and times
    times = [t for t in df['time'].unique()]
    id_choice = st.multiselect("Choose a cluster",
                               st.session_state['Idents'],
                               default=st.session_state['Idents'][1],
                               key="cluster_het")

    t_choice = st.selectbox('Pick the time', times, placeholder='ZT or CT...', key='time_het')

    # Select relevant data
    df = df[(df["Idents"].isin(id_choice)) & (df['time'] == t_choice)]
    df = df.select_dtypes(include=('float', 'int'))
    df = df - df.mean(axis=0)

    # Select colormap
    cmap = st.selectbox("Select colormap", COLORMAPS, key='cmap_het')

    # Plot
    fig = compile_heterogeneity(df, id_choice, t_choice, cmap)
    st.session_state['heatmap'] = fig


def make_matrixplots() -> None:
    """
    Design matrix plots (i.e: heatmaps) just like dotplots
    """

    # Backend data reordering
    adata = deepcopy(st.session_state['adata'])
    exps = ['LD', 'DD']

    st.write('## Select data to display')
    condition_choice = st.multiselect("Select condition", exps, default='DD', key="cond_mat")
    adata_plot = adata[adata.obs['condition'].isin(condition_choice)].copy()

    # Group by
    group = st.radio('Group by:', ['time', 'Idents'], key='groupby_mat')
    extra_str = ''
    if group == 'time':
        idents = st.multiselect("Select one or more clusters",
                                st.session_state['Idents'],
                                default=st.session_state['Idents'][0],
                                key="cluster_mat")
        adata_plot = adata_plot[adata_plot.obs['Idents'].isin(idents)].copy()
        extra_str = f" for {','.join(idents)}"

    # Choose gene subset
    var_names = st.multiselect("Determine gene order", adata_plot.var_names.unique(),
                               default=st.session_state['genes'],
                               key="vnames_mat")

    swap_axes = st.checkbox('Show groups in X axis', key='swap_axes_mat')

    # Select colormap
    cmap = st.selectbox("Select colormap", COLORMAPS, key='cmap_mat')

    # Plot!
    fig = compile_scplot(plot_type='matrixplot',
                         _adata_plot=adata_plot,
                         group=group,
                         condition_choice=condition_choice,
                         var_names=var_names,
                         extra_str=extra_str,
                         swap_axes=swap_axes,
                         cmap=cmap)

    st.session_state['matrixplot'] = fig


def make_tsne_plots() -> None:
    """
    Makes t-SNE plots for the selected data
    Returns:
    """
    # Select colormap
    cmap = st.selectbox("Select colormap", COLORMAPS, key='cmap_tsne')
    compile_tsne_plot(cmap)
    if st.session_state['tsne_plots'] is not None:
        for figure in st.session_state['tsne_plots']:
            st.pyplot(figure)


def make_umap_plots() -> None:
    """
    Makes UMAP plots for the selected data
    Returns:
    """
    # Select colormap
    cmap = st.selectbox("Select colormap", COLORMAPS, key='cmap_umap')
    compile_umap_plot(cmap)
    if st.session_state['umap_plots'] is not None:
        for figure in st.session_state['umap_plots']:
            st.pyplot(figure)


####################################################################################
#                                MAIN UI                                           #
####################################################################################

def main():
    color_theme_setup()

    # Page Title
    st.set_page_config(page_title="CircDrosView", page_icon="bar-chart")
    st.write("""
    # Interactive visualizer for *A transcriptomic taxonomy of* Drosophila *circadian neurons around the clock*
             """)

    # Initialization
    with st.spinner(text='Initializing variables...'):
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
            st.session_state['tsne_plots'] = None
            st.session_state['umap_plots'] = None
        fulladata = get_full_adata()
        all_genes = fulladata.var_names.tolist()

    # Gene selection
    st.write('## Step 1: select genes to analyze')
    with st.form(key="input_parameters"):
        choice = st.multiselect("Choose your genes", all_genes)
        submit = st.form_submit_button("Submit")

    if submit:
        fetch_data(choice)

    if not st.session_state['dataframe'].empty:
        # Download button for data
        df = st.session_state['dataframe']
        st.write("Download your data")
        st.download_button(label="Download filtered data as CSV",
                           data=df.to_csv().encode('utf-8'),
                           file_name='data.csv')

        # Tabs
        tab_dot, tab_point, tab_het, tab_mat, tab_tsne, tab_umap = st.tabs(['Dot plots',
                                                                            'Point plots',
                                                                            'Heterogeneity heatmap',
                                                                            'Matrix plot',
                                                                            't-SNE plots',
                                                                            'UMAP plots'])
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

        # Matrixplots
        with tab_mat:
            plt.close()
            make_matrixplots()
            if st.session_state['matrixplot'] is not None:
                st.pyplot(st.session_state['matrixplot'])

        # t-SNE
        with tab_tsne:
            plt.close()
            make_tsne_plots()
        # UMAP
        with tab_umap:
            plt.close()
            make_umap_plots()

    else:
        st.write("Please select data to fetch.")


if __name__ == '__main__':
    main()
