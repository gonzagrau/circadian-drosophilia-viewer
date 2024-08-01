import streamlit as st
import pandas as pd
import scanpy as sc
from dataset_handler import read_genexp_files


def main():
    ex_genes = ["DIP-gamma", "DIP-beta", "DIP-delta", "DIP-theta"]

    ######################
    # Page Title
    ######################

    st.write('''
    # Rosbash (2021) circadian *Drosophilia* data visualization tool
    ''')

    options = st.multiselect("Choose your genes", ex_genes)
    df = read_genexp_files(options)
    st.dataframe(df)

if __name__ == '__main__':
    main()