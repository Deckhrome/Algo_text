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
if os.path.exists('../GENOME_REPORTS/overview.txt') and os.path.isfile("../pickle/organism_df"):
    t = os.path.getmtime("../GENOME_REPORTS/overview.txt")


print("Suppression de l'arborescence actuelle...")
try: shutil.rmtree('../GENOME_REPORTS')
except: pass
print("terminé.")
os.mkdir("../GENOME_REPORTS")
os.chdir('../GENOME_REPORTS')
os.mkdir('IDS')

with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
    ftp.login()  # connecter au FTP

    ftp.cwd('genomes/GENOME_REPORTS')

    print(" - overview.txt")
    ftp.retrbinary('RETR overview.txt', open("overview.txt",'wb').write)

    ftp.cwd('IDS')  # changer de répertoire courant
    n = len(ftp.nlst())
    for i,filename in enumerate(ftp.nlst()):
        print(f" - {filename}")
        ftp.retrbinary('RETR '+ filename, open("IDS/" + filename, 'wb').write)



    # parse overview.txt
    organism_names = []
    organism_paths = []

    with open('../GENOME_REPORTS/overview.txt') as f:
        print("traitement overview.txt...")
        first_row = True
        count_rows = 1
        for row in f:
            count_rows += 1
            if first_row:
                first_row=False
                continue
            parsed_row = row.split('\t')

            try :
                organism = parsed_row[0].replace(' ','_').replace('/','_')
                kingdom = parsed_row[1].replace(' ','_').replace('/','_')
                group = parsed_row[2].replace(' ','_').replace('/','_')
                subgroup = parsed_row[3].replace(' ','_').replace('/','_')
                path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
                organism_names.append(parsed_row[0])
                organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/')
            except IndexError :
                print(f"IndexError sur {parsed_row}")
                pass
        print("terminé traitement overview.txt")

# parse fichier ids
ids_files = os.listdir('../GENOME_REPORTS/IDS/')
organism_names_ids = []
organism_paths_ids = []
organism_NC_ids = []
i = 0

for ids in ids_files:
    print(f"Traitement ids : {ids}")
    i += 1

    with open('../GENOME_REPORTS/IDS/' + ids) as f:
        n_line = sum(1 for _ in f)

    with open('../GENOME_REPORTS/IDS/' + ids) as f:
        for i,row in enumerate(f):
            # update l'affichage de la progression

            parsed_row = row.replace('\n', '').split('\t')
            if (parsed_row[1][0:2] != 'NC'):
                continue

            try:
                index = organism_names.index(parsed_row[5])

            except ValueError:
                # parsed_name = parsed_row[5].split(' ')[::-1]
                # try_name = parsed_row[5]
                # for word in parsed_name :
                #     try_name = try_name.replace(' '+word, '')
                #     try:
                #         index = organism_names.index(try_name)
                #         break
                #     except : pass
                continue # not found in overview.txt

            try:
                if parsed_row[1] not in organism_NC_ids[organism_names_ids.index(organism_names[index])]:
                    organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
            except ValueError:
                organism_names_ids.append(organism_names[index])
                organism_paths_ids.append(organism_paths[index])
                organism_NC_ids.append([parsed_row[1]])
                name = organism_names[index].replace(" ", "_")
                name = name.replace("[", "_")
                name = name.replace("]", "_")
                name = name.replace(":", "_")
                name = name.replace("/", "_")
                path = organism_paths[index] + name + "/"
                if not os.path.exists(path):
                    os.makedirs(path)



organism_df = pd.DataFrame({
            "name":organism_names_ids,
            "path":organism_paths_ids,
            "NC":organism_NC_ids})

# crée le fichier pickle sauvegardant le dataframe
if not os.path.exists("../pickle"):
    os.makedirs("../pickle")

with open("../pickle/organism_df", 'wb') as f:
    pickle.dump(organism_df, f)

df = pd.read_pickle("../pickle/organism_df")

df["nb_NC"] = df["NC"].apply(len)


homo = df[df['name'] == "Homo sapiens"]

print(homo)

NC_list = homo["NC"].iloc[0]
NC_list = list(dict.fromkeys(NC_list))

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
a=f(NC_list)