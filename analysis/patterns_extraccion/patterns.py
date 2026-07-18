#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
El objetivo de este codigo es extraer informacion de los archivos de salidad de Hmmer de cada gen, para luego crear una matriz de presencia/ausencia deg
genes de cada genoma.

Este script toma en cuenta que la salida de Hmmer tiene una estructura especifica, es decir, que esta tiene una columna de ID del genoma y
una columna final con la asignacion de la proteina. El script se enfoca en extraer la proteina que mas se repire en las primeras 30 lineas, y luego contar 
cuantas veces se encuentra esa proteina en cada genoma.

El output final es un archivo CSV con una matriz donde las filas son los genomas
"""
#Cargamos las librerias
import pandas as pd             # esta libreria para manipular dataframes
import re                       # esta ayuda a buscar patrones de texto
from pathlib import Path        # esto ayuda a manejar las rutas de los archivos
from collections import defaultdict # solo crea valores por defecto cuando no hay


# ruta de los ID de los genomas que quiero analizar 
RUTA_ARCHIVO_IDS = '/mnt/data/sur/users/mreyes/data/beggiatoas/results/checkm/storage/ID.txt'

# esto lee los IDs
def cargar_lista_ids(ruta):
  """
    Lo que hace esta funcion es leer el archivo ID.txt y extraer los numeros de ID de genomas que se encuentran en ese archivo.
     - El archivo ID.txt se espera que tenga una estructura donde cada linea contiene un ID de genoma.
    """
  ids = [] # lista vacia para guardar los IDs
  ruta_path = Path(ruta)

  # simplemente por seguridad y ver si el archivo ID se encuentra en la ruta.
  if not ruta_path.exists():
    print(f" No se encontró el archivo de IDs en: {ruta}")
    return []

  print(f"Leyendo IDs desde: {ruta} ...") #si se encontro

  with open(ruta_path, 'r', encoding='utf-8') as f: # se abre el archivo en modo lectura
    for linea in f:
      linea = linea.strip() # Quita espacios y saltos de línea al inicio/final
      if not linea: continue # Si la línea está vacía, se la salta

      partes = linea.split() # se divide la linea en partes, por ejemplo si la linea es "12345 some description", partes seria ["12345", "some", "description"]
      posible_id = partes[-1] # se asume que el ID del genoma es la última parte de la línea, por ejemplo "12345" en el ejemplo anterior

      #  si se cuentran numeros puede seguir
      if posible_id.isdigit():
        ids.append(int(posible_id))
      else:
        print(f"Ignorando línea no numérica: {linea}") #para ver si no se encontro una linea que no sea numerica

  return ids

# se ejecuta la funcion para cargar los IDs y se guarda en una variable global
LISTA_GENOMAS = cargar_lista_ids(RUTA_ARCHIVO_IDS)

# Si la lista está vacía,el programa se detiene
if not LISTA_GENOMAS:
  print("La lista de genomas está vacía")

#rutas de los archivos con los genes, ya que en HMMER cada gen se encuentra en un archivo diferente.
ARCHIVOS_GENES = {
  'murG': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/murG/resultados_murG.txt',
  'murC': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/murC/resultados_murC.txt',
  'murB': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/murB/resultados_murB.txt',
  'ddl': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/ddl/resultados_ddl.txt',
  'ftsZ': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/ftsZ/resultados_ftsZ.txt',
  'lpxC': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/lpxC/resultados_lpxC.txt',
  'ftsA': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/ftsA/resultados_ftsA.txt',
  'ftsQ': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/ftsQ/resultados_ftsQ.txt',
  'mreB': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/mreB/resultados_mreB.txt',
  'mreC': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/mreC/resultados_mreC.txt',
  'mrdA1': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/mrdA1/resultados_mrdA1.txt',
  'mrdA2': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/mrdA2/resultados_mrdA2.txt',
  'mltB': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/mltB/resultados_mltB.txt',
  'mreD1': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/mreD1/resultados_mreD1.txt',
  'rodZ': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/rodZ/resultados_rodZ.txt',
  'rodA': '/mnt/data/sur/users/mreyes/data/beggiatoas/results/aligments/rodA/resultados_rodA.txt'
}

# ruta de salida donde se guardaran los resultados finales.
RUTA_SALIDA = '/mnt/data/sur/users/mreyes/data/beggiatoas/results/resultado_Pattern.csv'

def extraer_proteina_principal(archivo, n_lineas=30):
  """
    Esta función lee las primeras n_lineas de un archivo HMMER y extrae la descripción de la proteína que más se repite.
      - Se asume que la descripción de la proteína se encuentra a partir de la columna 19 en adelante (index 18 en Python).
      - Se utiliza un contador para encontrar la proteína más común entre las primeras n_lineas.  
      - Si no se encuentra ninguna descripción válida, la función devuelve None.
    """
  descripciones = [] # lista para guardar las descripciones de las proteínas encontradas en las primeras n_lineas

  if not Path(archivo).exists():
    return None

  with open(archivo, 'r', encoding='utf-8', errors='ignore') as f: # se abre el archivo en modo lectura, con manejo de errores para caracteres no válidos
    for linea in f:
      if linea.startswith('#') or linea.strip() == '': # Si la línea es un comentario (comienza con '#') o está vacía, se la salta
        continue

      partes = linea.split() # lo que 
      if len(partes) < 19: #
        continue

      descripcion = ' '.join(partes[18:]) # se une la parte de la descripción de la proteína (a partir de la columna 19) en una sola cadena
      descripciones.append(descripcion) # se agrega la descripción a la lista de descripciones

      # se leen las primeras n lineas
      if len(descripciones) >= n_lineas:
        break

  if not descripciones:
    return None

  #counter servira para que se encuentre la moda de la proteina que mas se repite
  from collections import Counter
  proteina_principal = Counter(descripciones).most_common(1)[0][0]
  return proteina_principal


def parsear_hmmer_gen(archivo, lista_genomas):
  """
    El objetivo de esta funcion es leer un archivo de salida de HMMER para un gen específico, 
    identificar la proteína principal que se encuentra en las primeras líneas del archivo, y 
    luego contar cuántas veces esa proteína aparece en cada genoma de la lista proporcionada.

    """
  # se extrae la proteina principal del archivo, la proteina que mas se repite en las primeras n_lineas del archivo.
  proteina_principal = extraer_proteina_principal(archivo)

  if proteina_principal is None:
    print(f"No se encontro proteina principal en {Path(archivo).name}")
    return {}, {}

  print(f"  -> Proteina principal identificada: {proteina_principal[:40]}...") 

  conteo = defaultdict(int) # diccionario que guarda las proteínas
  otras_proteinas = defaultdict(int) # esto solo indica que se descarto,para luego verlo manual
  total_lineas = 0

# esta parte es la que hace el conteo de las proteínas, se lee el archivo linea por linea, se extrae el ID del genoma y la descripción de la proteína, 
# y se cuenta solo si la descripción coincide con la proteina principal.
  with open(archivo, 'r', encoding='utf-8', errors='ignore') as f: #
    for linea in f:
      if linea.startswith('#') or linea.strip() == '':
        continue

      partes = linea.split()
      if len(partes) < 19: continue

      # se extrae el Id del gen para saber que genoma es
      nombre_completo = partes[0]
      try:
        genome_str = nombre_completo.split('_')[0] 
        genome_id = int(genome_str)
      except ValueError:
        continue
      # solo se cuentan los genomas que estan en el ID.txt
      if genome_id not in lista_genomas:
        continue

      # se extrae la descripcion de cada linea
      descripcion = ' '.join(partes[18:])

      #se cuenta solo se la descripcion que coinciden con la proteina principal
      if descripcion == proteina_principal:
        conteo[genome_id] += 1
      else:
        otras_proteinas[descripcion] += 1

  # reporte final
  stats = {
    'proteina_principal': proteina_principal,
    'hits_principal': sum(conteo.values()),
    'otras_proteinas': dict(otras_proteinas),
    'total_lineas': total_lineas
  }

  return dict(conteo), stats


def construir_matriz_genes(archivos_dict, lista_genomas):
  """
  esta funcion construye una matriz de presencia/ausencia de genes para cada genoma, a partir de los archivos de salida de HMMER para cada gen.
  - Se inicializa una matriz con ceros, donde las filas son los genomas y las columnas son los genes.
  - Para cada gen, se procesa su archivo de salida de HMMER, se extrae la proteína principal y se cuenta cuántas veces aparece en cada genoma.
  - Se rellena la matriz con los conteos obtenidos para cada genoma y gen.
  - Se devuelve la matriz como un DataFrame de pandas y un diccionario con estadísticas sobre las proteínas encontradas.
    """
  # Inicializar matriz con ceros: Filas=Genomas, Columnas=Genes
  matriz = {genoma: {gen: 0 for gen in archivos_dict.keys()} 
    for genoma in lista_genomas}
  estadisticas = {}

  for gen, archivo in archivos_dict.items(): # se itera sobre cada gen y su archivo correspondiente
    print(f"\n Procesando gen: {gen}...") #alucionacion de claude que funciona, por que me dice como va

    if not Path(archivo).exists():
      print(f"Archivo no encontrado: {archivo}")
      continue

    # conteo del gen especifico
    conteo, stats = parsear_hmmer_gen(archivo, lista_genomas)

    # se rellena la matriz
    for genome_id, cantidad in conteo.items():
      if genome_id in matriz: 
        matriz[genome_id][gen] = cantidad

    estadisticas[gen] = stats
    print(f"   -> {len(conteo)} genomas tienen este gen.")

  # Convertir el diccionario a un dataframe 
  df = pd.DataFrame.from_dict(matriz, orient='index')
  df.index.name = 'genome_id'
  df.reset_index(inplace=True) # Saca el ID del índice y lo vuelve columna

  return df, estadisticas


def asignar_pattern1(row):
  """
  Esta funcion asigna un pattern a cada genoma basado en la presencia o ausencia de un conjunto específico de genes relacionados con la división celular.
    Pattern 1: Cluster completo de división celular.
    Regla: Deben estar TODOS los genes de la lista presentes (valor >= 1).
    """
  genes_pattern1 = ['murG', 'murC', 'murB', 'ddl', 'ftsZ', 
                    'lpxC', 'ftsA', 'ftsQ',  'mreB', 'mreC', 
                    'mrdA1', 'mrdA2', 'mltB', 'mreD1', 'rodZ']

  genes_existentes = [g for g in genes_pattern1 if g in row.index]

  # se verifica si todos tienen al menos una copia de cada gen
  todos_presentes = all(row[gen] >= 1 for gen in genes_existentes)

  if len(genes_existentes) < len(genes_pattern1):
    return '-'

  return '+' if todos_presentes else '-'


def asignar_pattern2(row):
  """
  Esta función asigna un pattern a cada genoma basado en la duplicación de ciertos genes del elongasoma y la ausencia de genes del divisoma.
    Pattern 2: Duplicación del 'elongasome' y ausencia de divisoma.
    Regla: (mreD1, mrdA1, rodZ >= 2 copias) Y (ftsQ, ftsA == 0 copias).
    """
  #se verifica la duplicacion de los genes:
  duplicados = (
    'mreD1' in row.index and row['mreD1'] >= 2 and
    'mrdA1' in row.index and row['mrdA1'] >= 2 and  
    'rodZ' in row.index and row['rodZ'] >= 2
  )

  # se verifica la ausencia de los genes:
  incompleto = (
    'ftsQ' in row.index and row['ftsQ'] == 0 and
    'ftsA' in row.index and row['ftsA'] == 0
  )

  # Deben cumplirse ambas condiciones para decir que estan duplicador
  return '+' if (duplicados and incompleto) else '-'

if __name__ == '__main__':
  print("=" * 70)
  print("EXTRACCION DE GENES Y ASIGNACION DE PATTERNS")
  print(f"IDs cargados para analizar: {len(LISTA_GENOMAS)}")
  print("=" * 70)

  if len(LISTA_GENOMAS) == 0:
    print(" No se puede continuar sin genomas")
    exit()

  # matriz de conteo de genes por genoma
  print("\n Construyendo matriz de genes...")
  df_genes, estadisticas = construir_matriz_genes(ARCHIVOS_GENES, LISTA_GENOMAS)

  # resumen en pantalla de cuantos genomas tienen cada gen, y cuantos tienen duplicaciones
  print("\n" + "=" * 70)
  print("RESUMEN DE GENES ENCONTRADOS")
  print("=" * 70)
  for gen in ARCHIVOS_GENES.keys():
    if gen in df_genes.columns:
      presentes = (df_genes[gen] > 0).sum()
      duplicados = (df_genes[gen] > 1).sum()
      print(f"{gen:8s}: {presentes:3d} genomas tienen el gen ({duplicados} tienen duplicaciones)")

  # asignar los patterns + / - dependiendo la situacion
  print("\n Asignando patterns (+/-)...")
  df_genes['pattern1'] = df_genes.apply(asignar_pattern1, axis=1)
  df_genes['pattern2'] = df_genes.apply(asignar_pattern2, axis=1)

  # se cuenta cuantos genes cumplieron y cuantos no cumplieron cada pattern
  print("\n" + "=" * 70)
  print("RESUMEN DE PATTERNS")
  print("=" * 70)
  print("\nPattern 1 (Cluster completo):")
  print(df_genes['pattern1'].value_counts())

  print("\nPattern 2 (Elongasoma duplicado + Divisoma incompleto):")
  print(df_genes['pattern2'].value_counts())

  # Archivo sencillo para R, para input de la grafica, solo con el ID del genoma y los patterns asignados
  df_salida = df_genes[['genome_id', 'pattern1', 'pattern2']].sort_values('genome_id')
  df_salida.to_csv(RUTA_SALIDA, index=False, header=False)
  print(f"\n Archivo final para gráfica guardado en: {RUTA_SALIDA}")

  # Archivo detallado con la matriz completa de genes, para que se pueda revisar manualmente si se quiere, con toda la informacion de los conteos de cada gen.
  ruta_detallada = RUTA_SALIDA.replace('.csv', '_detallado.csv')
  df_genes.to_csv(ruta_detallada, index=False)
  print(f" Archivo detallado (matriz completa) guardado en: {ruta_detallada}")

  # Archivo de validacion con un resumen de las proteinas principales encontradas en cada gen, y cuantos hits se contaron para cada una, para que se pueda revisar manualmente si se quiere.
  ruta_validacion = RUTA_SALIDA.replace('.csv', '_validacion.txt')
  with open(ruta_validacion, 'w', encoding='utf-8') as f:
    f.write("REPORTE DE VALIDACION\n")
    f.write(f"Total Genomas analizados: {len(LISTA_GENOMAS)}\n\n")
    for gen, stats in estadisticas.items():
      f.write(f"GEN: {gen}\n")
      f.write(f"  Proteina principal detectada: {stats.get('proteina_principal', 'N/A')}\n")
      f.write(f"  Total Hits contados: {stats.get('hits_principal', 0)}\n")
      f.write("-" * 30 + "\n")

  print("\nProceso finalizado con éxito.") # indica que el proceso se ha completado sin errores
