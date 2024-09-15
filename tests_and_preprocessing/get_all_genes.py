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
    LD_path = r'data_subsets/GSM4768102_ZT22_20190309_AR01.csv'
    LD_genes = list_indeces(LD_path)
    LD_genes.sort()
    DD_path = r'data_subsets/GSM4768057_CT18_20190710_AR18.csv'
    DD_genes = list_indeces(DD_path)
    DD_genes.sort()

    all_genes = set(LD_genes).intersection(set(DD_genes))
    all_genes = list(all_genes)
    all_genes.sort()
    print(len(all_genes))

    with open('../LD_genes.txt', 'w') as file:
        for gene in LD_genes:
            file.write(f"{gene}\n")

    with open('../DD_genes.txt', 'w') as file:
        for gene in DD_genes:
            file.write(f"{gene}\n")

    with open('../all_genes.txt', 'w') as file:
        for gene in all_genes:
            file.write(f"{gene}\n")


if __name__ == '__main__':
    main()