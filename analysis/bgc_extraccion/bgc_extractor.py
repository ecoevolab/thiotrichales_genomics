'''
El objetivo de este script es extraer el numero de pares de bases que ocupan los BGCs predichos
por antismash, asi como el producto que se predice para cada uno de ellos. A traves de un proceso de iteracion donde
se leen los archivos .gbk de cada region predicha, se extraen los datos necesarios y se guardan en un archivo .csv con el siguiente formato:
genome_id, bgc_id, start_bp, end_bp, length_bp, product

El script se ejecuta desde la terminal con el siguiente comando:
python bgc_extractor.py --folder /ruta/a/la/carpeta/antismash --out /ruta/al/archivo/salida.csv
'''

#Cargar librerias 
import argparse
import csv
import re
from pathlib import Path


def parse_region_gbk(gbk_file): 
    '''
    la finalidad de esta funcion es tomar un archivo .gdk de una region predicha por antismash y extraer la siguiente informacion:
    - start_bp: la posicion de inicio del BGC en la secuencia genómica
    - end_bp: la posicion de fin del BGC en la secuencia genómica
    - length_bp: la longitud del BGC en pares de bases (end_bp - start_bp + 1)
    - product: el producto que se predice para el BGC, si esta disponible, o "NA" si no se encuentra la informacion del producto en el archivo .gbk
    '''
    start = end = product = "NA" #Se inicia todo en NA para manejar casos donde no se encuentre la informacion en el archivo .gbk

    with open(gbk_file) as f: 
        for line in f:
            # región principal
            m = re.match(r'\s*region\s+(\d+)\.\.(\d+)', line) # se busca la linea que contiene la informacion de la region principal del BGC, que tiene el formato "region start..end"
            if m:
                start, end = map(int, m.groups()) #caso exitoso, se extrae el inicio y el fin del BGC y se convierten a enteros

            # producto
            if '/product=' in line:
                product = line.split('=')[1].replace('"', '').strip() # se busca la linea que contiene la informacion del producto y se extrae el nombre del producto, eliminando comillas y espacios en blanco

    if start != "NA" and end != "NA": #si se encuentra la informacion de inicio y fin del BGC, se calcula la longitud del BGC en pares de bases
        length = end - start + 1 
        return start, end, length, product

    return None #caso donde no se encuentra la informacion de inicio y fin del BGC


def process_antismash(parent_folder, output_csv): 
    ''' Esta funcion se encarga de procesar la carpeta que contiene los resultados de antismash, 
    iterando sobre cada subcarpeta (cada una correspondiente a un genoma) y 
    extrayendo la informacion de cada region predicha por antismash utilizando la funcion parse_region_gbk. 
        '''

    with open(output_csv, "w", newline="") as out: 
        writer = csv.writer(out) #se crea un objeto writer para escribir en el archivo .csv de salida
        writer.writerow([ #formato del encabezado del archivo .csv de salida
            "genome_id",
            "bgc_id",
            "start_bp",
            "end_bp",
            "length_bp",
            "product"
        ])

        for genome_dir in Path(parent_folder).iterdir(): #se itera sobre cada subcarpeta en la carpeta principal
            if not genome_dir.is_dir(): #caso donde no se es una carpeta
                continue

            genome_id = genome_dir.name #se obtiene el nombre de la carpeta y se asigna como genome_id
            bgc_id = 1 #se inicia el contador de bgc_id en 1 para cada genoma

            region_files = sorted(genome_dir.glob("*.region*.gbk")) #se buscan los archivos .gbk que corresponden a las regiones predichas por antismash

            for gbk in region_files: #se itera sobre cada archivo .gbk encontrado
                data = parse_region_gbk(gbk) # se llama a la funcion parse_region_gbk para extraer la informacion del archivo .gbk
                if data:
                    writer.writerow([ #se escribe una fila en el archivo .csv de salida
                        genome_id,
                        bgc_id,
                        *data
                    ])
                    bgc_id += 1


if __name__ == "__main__": 
    parser = argparse.ArgumentParser() 
    parser.add_argument("--folder", required=True) 
    parser.add_argument("--out", required=True) 
    args = parser.parse_args() 

    process_antismash(args.folder, args.out) #se llama a la funcion process_antismash con los argumentos proporcionados por el usuario a través de la terminal

