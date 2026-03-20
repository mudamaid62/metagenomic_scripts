import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches # Para la leyenda manual
import sys
import warnings

# --- Configuración de Archivos y Parámetros ---
DATA_FILE = "metadata_con_umap.tsv"
OUTPUT_PDF = "umap_scatter_pident_centroids.pdf"

# Parámetros gráficos
S_VAL = 10
LW_VAL = 0.1
ALPHA_BASE = 0.2     # Transparencia para la capa viridis (fondo)
ALPHA_TOP = 0.7      # Transparencia para los centroides (resaltado)
FIG_SIZE = (12, 10)

# --- 1. Carga de Datos (con Pandas) ---
print(f"Cargando datos desde: {DATA_FILE}")
try:
    df = pd.read_csv(DATA_FILE, sep='\t')
    
    # Verificar columnas necesarias
    required_cols = ["UMAP_1", "UMAP_2", "Single", "Centroid", "Pident_modifier"]
    if not all(col in df.columns for col in required_cols):
        print("Error: Faltan columnas. Se requieren:")
        print(required_cols)
        print(f"Columnas disponibles: {df.columns}")
        sys.exit(1)
        
    print(f"-> Datos cargados: {df.shape[0]} filas.")

except FileNotFoundError:
    print(f"Error: No se encontró el archivo '{DATA_FILE}'.")
    sys.exit(1)
except Exception as e:
    print(f"Error al cargar el archivo: {e}")
    sys.exit(1)

# --- 2. Filtrado y Preparación de Datos ---
print("Filtrando datos...")

# Paso 1: Excluir 'Single == "Yes"'
df_filtrado = df[df['Single'] == 'No'].copy()
rows_dropped = len(df) - len(df_filtrado)
print(f"-> Se excluyeron {rows_dropped} filas ('Single == Yes').")
print(f"-> Quedan {len(df_filtrado)} filas para graficar.")

# Paso 2: Crear subconjuntos para las capas
df_base = df_filtrado[df_filtrado['Centroid'] == 'No']
df_top = df_filtrado[df_filtrado['Centroid'] == 'Yes']

print(f"-> Capa Base ('Centroid == No'): {len(df_base)} puntos (Viridis)")
print(f"-> Capa Superior ('Centroid == Yes'): {len(df_top)} puntos (Rojo)")


# --- 3. Creación del Gráfico (Capas Múltiples) ---
print(f"Generando Gráfico -> {OUTPUT_PDF} ...")
try:
    fig = plt.figure(figsize=FIG_SIZE)
    #sns.set_theme(style="ticks")
    ax = fig.add_subplot(111)
    
    # --- Capa 1: Base (Pident_modifier - Viridis) ---
    print("Dibujando capa base (Pident)...")
    g = sns.scatterplot(
        data=df_base,
        x="UMAP_1",
        y="UMAP_2",
        hue="Pident_modifier", # Gradiente según esta columna
        palette="viridis",     # Paleta Viridis
        s=S_VAL,
        linewidth=0,
        alpha=ALPHA_BASE,      # Transparencia base
        legend="auto"          # Crear la barra de color
    )
    
    # --- Capa 2: Superior (Centroides - Rojo) ---
    if not df_top.empty:
        print("Dibujando capa superior (Centroides)...")
        sns.scatterplot(
            data=df_top,
            x="UMAP_1",
            y="UMAP_2",
            color="red",         # Color rojo sólido
            s=S_VAL,
            linewidth=0,
            alpha=ALPHA_TOP,     # Más opaco para resaltar
            legend=False,        # No crear leyenda para esta capa
            ax=g # IMPORTANTE: Plotea en el mismo eje (g)
        )
        
    #plt.title("UMAP: Pident (Viridis) y Centroides (Rojo)", fontsize=18)
    #plt.xlabel("Componente UMAP 1", fontsize=12)
    #plt.ylabel("Componente UMAP 2", fontsize=12)
    
    # Añadir títulos y etiquetas
    #plt.title('pident', fontsize=16)
    plt.xlabel('UMAP_1', fontsize=12)
    plt.ylabel('UMAP_2', fontsize=12)

    # Configurar la cuadrícula: líneas grises cada 2 unidades
    ax.grid(True, which='major', axis='both', color='gray', linewidth=0.5, alpha=0.7)
    ax.set_xticks(range(int(ax.get_xlim()[0]) // 2 * 2, int(ax.get_xlim()[1]) + 2, 2))  # Múltiplos de 2 en X
    ax.set_yticks(range(int(ax.get_ylim()[0]) // 2 * 2, int(ax.get_ylim()[1]) + 2, 2))  # Múltiplos de 2 en Y

    # Opcional: mejorar apariencia de los ticks
    ax.tick_params(axis='both', which='major', labelsize=10)


    # --- Manejo de Leyendas ---
    
    # 1. Mover la barra de color (del 'hue')
    legend_colorbar = g.get_legend()
    if legend_colorbar:
        legend_colorbar.set_title("Pident_modifier")
        sns.move_legend(g, "upper left", bbox_to_anchor=(1.02, 1), frameon=False)

    # 2. Crear una leyenda manual para el punto rojo
    red_patch = mpatches.Patch(
        color='red', 
        alpha=ALPHA_TOP, 
        label='Centroid (Yes)'
    )
    
    # Añadimos la segunda leyenda manualmente,
    # en una posición diferente para que no choque
    ax = plt.gca() # Obtener el eje actual
    legend2 = ax.legend(
        handles=[red_patch],
        title="Capas",
        loc="lower left",
        bbox_to_anchor=(1.02, 0.0), # Abajo a la derecha
        frameon=False
    )
    # Volvemos a añadir la primera leyenda (la barra de color)
    # para que Matplotlib no la oculte.
    if legend_colorbar:
        ax.add_artist(legend_colorbar)

    sns.despine()
    
    # Guardar como PDF
    plt.savefig(OUTPUT_PDF, format="pdf", bbox_inches='tight')
    plt.close() # Cerrar la figura
    
    print(f"-> Gráfico guardado exitosamente en {OUTPUT_PDF}")

except Exception as e:
    print(f"Error al generar el gráfico: {e}")

print("\n--- ¡Proceso Completado! ---")
