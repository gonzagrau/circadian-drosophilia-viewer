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
    LD_path = r'..\..\dataset\GSM4768068_ZT02_20190309_AR08.csv'
    LD_genes = list_indeces(LD_path)
    LD_genes.sort()
    DD_path = r'..\..\dataset\GSM4768021_CT02_20190528_AR06.csv'
    DD_genes = list_indeces(DD_path)
    DD_genes.sort()

    all_genes = set(LD_genes).union(set(DD_genes))
    all_genes = list(all_genes)
    all_genes.sort()
    print(len(all_genes))

    with open('../database connection/LD_genes.txt', 'w') as file:
        for gene in LD_genes:
            file.write(f"{gene}\n")

    with open('../database connection/DD_genes.txt', 'w') as file:
        for gene in DD_genes:
            file.write(f"{gene}\n")

    with open('../all_genes.txt', 'w') as file:
        for gene in all_genes:
            file.write(f"{gene}\n")


if __name__ == '__main__':
    main()