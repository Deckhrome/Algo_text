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
print(max(df["nb_NC"]))

homo = df[df['name'] == "Culex"]

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