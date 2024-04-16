#!/usr/bin/env python3
from Bio import Entrez
from lxml import etree
import pandas as pd
import os
import shutil
from Bio import SeqIO
import pickle
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import random
import string
import signal
from ftplib import FTP
import re
import time
import os.path
import datetime
from threading import Thread
import functools

save_pickle = False


def timeout(timeout):
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [Exception('function [%s] timeout [%s seconds] exceeded!' % (func.__name__, timeout))]
            def newFunc():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e
            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(timeout)
            except Exception as je:
                print ('error starting thread')
                raise je
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco

def reset_tree(root, progressbar, lbl, progress = None, window = None):
    """
    Recrée l'arborescence stockée localement
    
    Si le fichier pickle existe est qu'il date d'aujourd'hui, alors on n'a pas besoin de le reset
    Sinon on essaye de le reset
    puis on recrée les dossiers
    et on recharge la liste des IDs de Genbank
    on convertit tout
    Avec les différentes classe, on fait l'arboressence
    Ensuite avec les différents noms on peut déjà préconstruire les dossiers pour chaque famille
    """
    

    if os.path.exists('../GENOME_REPORTS/overview.txt') and os.path.isfile("../pickle/organism_df"):
        t = os.path.getmtime("../GENOME_REPORTS/overview.txt")
        pred_time = datetime.datetime.fromtimestamp(t)
        now = datetime.datetime.now()

        if pred_time.month == now.month and pred_time.day == now.day and pred_time.year == now.year : # on update pas l'arbre
            root.destroy()
            return

    print("Suppression de l'arborescence actuelle...")
    lbl["text"] = "Suppression de l'arborescence actuelle..."
    try: shutil.rmtree('../GENOME_REPORTS')
    except: pass
    print("terminé.")
    os.mkdir("../GENOME_REPORTS")
    os.chdir('../GENOME_REPORTS')
    os.mkdir('IDS')

    print("Récupération des fichiers IDS...")
    lbl["text"] = f"Récupération de overview.txt"
    root.update()
    with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
        ftp.login()  # connecter au FTP

        ftp.cwd('genomes/GENOME_REPORTS')
        
        print(" - overview.txt")
        ftp.retrbinary('RETR overview.txt', open("overview.txt",'wb').write) 

        ftp.cwd('IDS')  # changer de répertoire courant
        n = len(ftp.nlst())
        lbl["text"] = f"Récupération des fichiers IDS (0/{n})"
        for i,filename in enumerate(ftp.nlst()):
            lbl["text"] = f"Récupération des fichiers IDS ({i+1}/{n})"
            progressbar["value"] = i/n * 100
            root.update()
            print(f" - {filename}")
            ftp.retrbinary('RETR '+ filename, open("IDS/" + filename, 'wb').write)

    # parse overview.txt
    organism_names = []
    organism_paths = []

    os.chdir('../script')
    with open('../GENOME_REPORTS/overview.txt') as f:
        print("traitement overview.txt...")
        lbl["text"] = "Traitement de 'overview.txt'"
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
        root.update()
        print(f"Traitement ids : {ids}")
        lbl["text"] = f"Traitement de '{ids}'"
        i += 1
        progressbar['value'] = 0
        root.update_idletasks()

        with open('../GENOME_REPORTS/IDS/' + ids) as f:
            n_line = sum(1 for _ in f)

        with open('../GENOME_REPORTS/IDS/' + ids) as f:
            for i,row in enumerate(f):
                # update l'affichage de la progression
                if i%(max(1,n_line//100)) == 0:
                    progressbar['value'] = i/n_line*100
                    root.update_idletasks()
                    root.update()

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

    root.destroy()



def check_inf_sup(inf,sup):
    """Fonction simple pour savoir qui est inferieur a qui

    Args:
        inf: value
        sup: value

    Returns:
        true if inf is inferior or equal to sup
        false if not
    """
    if(inf <= sup):
        return True
    else:
        return False
