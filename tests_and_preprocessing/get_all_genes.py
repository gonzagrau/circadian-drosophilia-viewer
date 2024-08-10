import pandas as pd


def list_indeces(filepath: str) -> list[str]:
    """
    Reads a csv file and returns a list of all genes names.
    Args:
        filepath: csv file path

    Returns:
        list of all genes names.
    """
    df = pd.read_csv(filepath, index_col=0, keep_default_na=False)
    index_df = df.index.to_list()
    return index_df


def main():
    zt_path = r'..\..\dataset\GSM4768068_ZT02_20190309_AR08.csv'
    zt_genes = list_indeces(zt_path)
    ct_path = r'..\..\dataset\GSM4768021_CT02_20190528_AR06.csv'
    ct_genes = list_indeces(ct_path)

    all_genes = set(zt_genes).union(set(ct_genes))
    all_genes = list(all_genes)
    all_genes.sort()
    print(len(all_genes))

    with open('../all_genes.txt', 'w') as file:
        for gene in all_genes:
            file.write(f"{gene}\n")


if __name__ == '__main__':
    main()