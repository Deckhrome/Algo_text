from ftplib import FTP
import os
import shutil
import pandas as pd
import pickle
import random
import string
import time
import requests
from Bio import Entrez
from Bio import SeqIO

df = pd.read_pickle("../pickle/organism_df")

df["nb_NC"] = df["NC"].apply(len)


homo = df[df['name'] == "Homo sapiens"]

NC_list = homo["NC"].iloc[0]
NC_list = list(dict.fromkeys(NC_list))

print(NC_list[:10])

def f(NC):

    Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20))+'@'+''.join(random.choice(string.ascii_lowercase) for i in range(20))+ '.com'
    try :
        print(f"  request for {NC} text")
        handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
    except Exception as e:
        print(f"1 for {NC} text",str(e))
        if str(e) == "HTTP Error 429: Too Many Requests":
            iter = 0
            iter_max = 10
            while True:
                time.sleep(10) # on se calme un peu
                try :
                    print(f"  request for {NC} text (rerun 1)")
                    handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
                    break
                except Exception as e:
                    iter += 1
                    print(f"1.{iter} for {NC} text",str(e))
                    if str(e) == "HTTP Error 429: Too Many Requests":
                        continue
                    elif iter > iter_max:
                        return
        else:
            return
    try:
        # record_fasta = timeout(timeout=TIMEOUT_MAX)(SeqIO.read)(handle_fasta, "fasta")
        record_fasta = SeqIO.read(handle_fasta, "fasta")
    except Exception as e:
        print(f"2 for {NC} text",e)
        return
    handle_fasta.close()

def process_NC(NC, path, name, selected_region):
    """
    Traite un identifiant NC en récupérant les données correspondantes depuis NCBI.

    Arguments :
    NC : Identifiant NC à traiter (string)
    path : Chemin de destination pour enregistrer les fichiers (string)
    name : Nom associé à l'identifiant NC (string)
    selected_region : Région spécifique à extraire (string)

    Renvoie :
    tuple : (résultat, code)
        - résultat : description du succès ou de l'échec du traitement (string)
        - code : couleur d'affichage du résultat (string)
    """
    Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '@' + ''.join(random.choice(string.ascii_lowercase) for i in range(20)) + '.com'

    try:
        print(f"Requesting data for {NC}...")
        handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
        record_fasta = SeqIO.read(handle_fasta, "fasta")
        handle_fasta.close()

        NC_filename = f"{selected_region}_{name}_{NC}.txt"

        with open(f"{path}/{name}/{NC_filename}", 'w') as out:
            for feature in record_fasta.features:
                if feature.type == selected_region:
                    out.write(f">{feature.id}\n{feature.location}\n{feature.qualifiers['translation'][0]}\n")

        return f"Successfully processed {NC}", "green"

    except Exception as e:
        print(f"Error processing {NC}: {str(e)}")
        return f"Error processing {NC}", "red"


# Exemple d'utilisation pour traiter un identifiant NC spécifique
NC_to_process = "NC_000001"
destination_path = "../test_save_data"
organism_name = "Homo_sapiens"
region_of_interest = "CDS"

result, color = process_NC(NC_to_process, destination_path, organism_name, region_of_interest)
print(result)
