import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    DATA_PATH = r'../dataset/GSM4768020_CT02_20190528_AR05.csv'
    ANNOT_PATH = r'../clock_neurons_annotation.csv'
    genexp_df = pd.read_csv(DATA_PATH, index_col=0).T
    print(genexp_df.head(10))
    annot_df = pd.read_csv(ANNOT_PATH, index_col=0)
    genexp_df['cluster'] = annot_df['Idents']


if __name__ == '__main__':
    main()