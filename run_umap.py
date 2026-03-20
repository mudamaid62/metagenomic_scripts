import pandas as pd
import numpy as np
import umap.umap as umap

#Carga tus datos en un DataFrame (reemplaza 'tus_datos.csv' con tu archivo)
df = pd.read_csv('iteration_4_go_matrix',sep="\t", index_col=0, header=0)

#Crea el objeto UMAP
reductor = umap.UMAP(n_components=2, random_state=42, n_neighbors=100, min_dist=0.01, metric='cosine')

#Aplica el UMAP al conjunto de datos
umap_datos = reductor.fit_transform(df.iloc[:, 1:])
#print(umap_datos)

#Añade los valores de UMAP a la dataframe
df_umap = pd.DataFrame(umap_datos, columns=['UMAP_1', 'UMAP_2'])
df_umap.to_csv('umap_output_v29102025.csv', header=False, index=False)
