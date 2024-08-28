import os
import pandas as pd
from dataset_handler import read_genexp_files
import anndata as ad

ANNOT_PATH = "C:/Users/Tassara/Documents/Doctorado/circadian-drosophilia-viewer-main/neuron_annotations.csv"
data_path: str = 'C:/Users/Tassara/Documents/Doctorado/circadian-drosophilia-viewer-main/dataset_LD/'
# Cambiar el nombre del nuevo dataset generado
anndata_PATH = "C:/Users/Tassara/Documents/Doctorado/circadian-drosophilia-viewer-main/dataset_LD2.h5ad"

annot_df = pd.read_csv(ANNOT_PATH, index_col=0)

genes = None
genexp_df = None
print('Going through files...')
for filename in os.listdir(data_path):
    new_df = pd.read_csv(os.path.join(data_path, filename), index_col=0) #, keep_default_na=False)
    new_df = new_df.rename(lambda x: x[1:], axis='columns').T  # removes leading 'x' char in idx strings
    if genes:
        try:
            new_df = new_df[genes]
        except KeyError:
            continue
    if genexp_df is None:
        genexp_df = new_df
        continue
    genexp_df = pd.concat([genexp_df, new_df])


full_df = genexp_df.join(annot_df, how='inner')

#print(full_df.shape)
#print(full_df.head(10))

#annot_cols = ['experiment', 'Repeats', 'condition', 'date', 'time', 'Idents']
#assert all([field in full_df.columns for field in annot_cols])
#obs = full_df[annot_cols]
#df = full_df.drop(annot_cols, axis='columns')

expression_matrix = full_df.iloc[:, :-6]  # Excluyendo las últimas 6 columnas
metadata = full_df.iloc[:, -6:]  # Solo las últimas 6 columnas son metadatos

adata = ad.AnnData(X=expression_matrix, obs=metadata)
adata.var['gene_names'] = expression_matrix.columns

# Guardar el objeto AnnData en un archivo .h5ad
adata.write(anndata_PATH)
print(adata)
print(adata.obs.head())
print(adata.var['gene_names'].head())