import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    DATA_PATH = r'../dataset/GSM4768020_CT02_20190528_AR05.csv'
    df = pd.read_csv(DATA_PATH)
    print(df.head(10))
    ax = sns.heatmap(data=df.head(10)[df.columns[1:]])
    ax.set_title('Ejemplo de lectura')
    plt.show()


if __name__ == '__main__':
    main()