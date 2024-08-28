import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata as ad

mpl.use('agg')
sns.set_palette('deep')

adata = sc.read_h5ad('dataset_LD.h5ad')
# Sacar genes
patterns_to_exclude = ['mt', 'rpls', 'rRNA', 'tRNA', 'ERCC', 'EGFP']
gene_names = adata.var['gene_names'].astype(str)
filtered_genes = ~gene_names.str.contains('|'.join(patterns_to_exclude), regex=True)

# Aplicar el filtro al objeto AnnData
adata_filtered = adata[:, filtered_genes].copy()
#print(adata_filtered.var['gene_names'])

# Normalizacion por celula (para que cada celula tenga el mismo total de cuentas)
# sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.normalize_total(adata_filtered, target_sum=1e4, exclude_highly_expressed=True)

# Aplicar logaritmo a los datos normalizados
#sc.pp.log1p(adata)

adata_filtered.obs['Idents'] = adata_filtered.obs['Idents'].astype('category')
print(adata_filtered.obs['Idents'])
# Seleccionar los genes de interes
#genes_of_interest = ['DIP-gamma', 'DIP-beta', 'DIP-delta', 'DIP-theta', 'dpr8']

# Crear el dotplot
mpl.use('TkAgg')  # Cambia a un backend interactivo, como TkAgg o Qt5Agg
#sc.pl.dotplot(adata, var_names=genes_of_interest, groupby='Idents', standard_scale='var',)

# Filtrar las celulas que pertenecen al cluster 'sLNVs'
adata_slNvs = adata_filtered[adata_filtered.obs['Idents'] == '2:s_LNv', :]
#adata_llNvs = adata_filtered[adata_filtered.obs['Idents'] == '25:l_LNv', :]

# Obtener los datos de expresion del gen
expression = adata_slNvs[:, adata_slNvs.var['gene_names'] == 'per'].X

# Crear un DataFrame para agrupar por 'time'
df = pd.DataFrame(expression.toarray(), columns=['per'], index=adata_slNvs.obs.index)
df['time'] = adata_slNvs.obs['time'].values

# Ordenar los puntos del tiempo
time_order = ['ZT02', 'ZT06', 'ZT10', 'ZT14', 'ZT18', 'ZT22']

# Crear el grafico usando Seaborn pointplot
plt.figure(figsize=(8, 5))
sns.pointplot(data=df, x='time', y='per', order=time_order, estimator='mean', errorbar='se')
plt.xlabel('Time (ZT)')
plt.ylabel('Mean Expression of per')
plt.title('Mean Expression in sLNVs Cluster Over Time')
plt.grid(True)
plt.show()
