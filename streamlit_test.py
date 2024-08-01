import streamlit as st
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_agg import RendererAgg
import pandas as pd
import scanpy as sc
from dataset_handler import read_genexp_files, df_to_anndata
from typing import List


mpl.use('TKAgg')


def fetch_data(choice: List[str]) -> None:
    """
    Reads data from csv files and stores it in session state
    Args:
        choice (List[str]): choice of genes to select from dataset
    """
    print(choice)
    with st.spinner(text='Fetching data...'):
        df = read_genexp_files(choice)
        st.success('Done')
    st.session_state['data'] = df



def make_dotplots(df: pd.DataFrame) -> None:
    """
    Interactive dotplot builder

    """
    adata = df_to_anndata(df)
    idents = [ident for ident in adata.obs['Idents'].unique()]
    idents.sort(key= lambda x : int(x.split(':')[0]))

    st.write('## Step 2: design your dotplots')
    group = st.radio('Group by:', ['exp_time','Idents'])
    adata_dotplot = adata
    if group == 'exp_time':
        idents = st.multiselect("Select a cluster", idents, default=idents[0])
        adata_dotplot = df_to_anndata(df[df['Idents'].isin(idents)])

    dotplot = sc.pl.DotPlot(adata_dotplot,
                            var_names=adata.var_names,
                            groupby=group,
                            standard_scale='var',
                            vmin=-1,
                            vmax=1,
                            var_group_rotation=0.,
                            edgecolors=None,
                            mean_only_expressed=True,
                            title=f"Gene expression by {group}",
                            cmap='Reds',
                            linewidth=0.)
    if group == 'exp_time':
        dotplot.swap_axes()
    dotplot.make_figure()
    st.pyplot(dotplot.fig)



def main():

    # Initialization
    st.set_page_config(page_title = "DRosviewer", page_icon="bar-chart")

    with open('all_genes.csv', 'r') as f:
        genes = f.read().splitlines()[1:]

    if 'data' not in st.session_state:
        st.session_state['data'] = pd.DataFrame()

    # Page Title
    TITLE = "# Interactive visualizer for *A transcriptomic taxonomy of* Drosophila *circadian neurons around the clock*"
    st.write(TITLE)


    # Gene selection
    st.write('## Step 1: select genes to analyze')

    with st.form(key="input_parameters"):
        choice = st.multiselect("Choose your genes", genes)
        submit = st.form_submit_button("Submit")

    if submit:
        fetch_data(choice)
        df = st.session_state['data']
        st.dataframe(df)
        make_dotplots(df)

    if not st.session_state['data'].empty:
        st.write("Fetched Data:")
        st.dataframe(st.session_state['data'])
        df = st.session_state['data']
        st.dataframe(df)
        make_dotplots(df)
    else:
        st.write("Please select data to fetch.")


if __name__ == '__main__':
    main()