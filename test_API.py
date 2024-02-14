# Importation des modules nécessaires de BioPython
from Bio import Entrez, SeqIO

# Fournir votre adresse e-mail pour l'utilisation de Entrez
Entrez.email = "jean.raphael2810@gmail.com"  # Remplacez par votre adresse e-mail réelle

# Liste des identifiants GenBank à parcourir (exemple fictif)
genbank_ids = ["NC_005816", "NC_007346", "NC_009899"]  # Remplacez ces valeurs par vos identifiants réels

# Boucle sur chaque identifiant GenBank dans la liste
for genbank_id in genbank_ids:
    print(f"Traitement de l'identifiant GenBank : {genbank_id}")

    # Utilisation de Entrez.efetch pour télécharger les données de séquence
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=genbank_id) as handle:
        seq_record = SeqIO.read(handle, "gb")  # Lecture de la séquence au format GenBank

        # Parcourir les features et extraire les CDS
        for feature in seq_record.features:
            if feature.type == "CDS":
                # Afficher l'identifiant de la CDS et sa séquence
                cds_id = feature.qualifiers.get('protein_id', [''])[0]  # Utiliser protein_id comme identifiant de CDS, si disponible
                cds_sequence = feature.extract(seq_record.seq)  # Extraire la séquence CDS à partir de la séquence complète
                print(f"CDS ID: {cds_id}")
                #print(f"Séquence CDS: {cds_sequence}\n")
                # Vous pouvez également afficher d'autres informations disponibles dans les qualifiers de la feature

    print(f"Fin du traitement pour l'identifiant GenBank : {genbank_id}\n---\n")

# Notez que ce script peut prendre un certain temps pour s'exécuter en fonction du nombre et de la taille des séquences à télécharger.
