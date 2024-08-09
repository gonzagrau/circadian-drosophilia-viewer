import pandas as pd


def main():
    PATH = r'..\..\dataset\GSM4768021_CT02_20190528_AR06.csv'
    df = pd.read_csv(PATH, index_col=0)
    index_df = df.index.to_frame()
    index_df.to_csv('all_genes.csv', index=False)


if __name__ == '__main__':
    main()