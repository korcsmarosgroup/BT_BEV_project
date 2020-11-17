# Predictinig microbe-host PPIs based on domain-motif interactions from ELM
# Input files:
    # Fasta sequences of human proteins which wants to be connected to microbes
    # motif details from ELM
    # domain - motif interactions from ELM
    # List of bacterial proteins and their domains

# Output file = all of the possible connections without quality filter
# Authors: Lejla Gul, David Fazekas

from pyfasta import Fasta
import re

def rename(fasta_key):
    fasta_key = fasta_key.split("|")
    fasta_key = fasta_key[1]
    return fasta_key


# fasta processing -> human.keys() print the keys, human[key_name] print the sequence
human = Fasta('/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/BT_OMV/resources/healthy/macrophage_healthy.fasta')


#elm identifier(key) - regex(value) dictionary
with open ("elm_classes_2020.tsv", "r") as motif_table:
    motif_table.readline()
    elm_regex = {}
    for line in motif_table:
        line = line.strip().split("\t")
        elm_regex[line[1]] = line[4]


#motif(key) - domain(value list) dictionary
with open ("elm_interaction_domains_2020.tsv", "r") as motif_domain_table:
    motif_domain_table.readline()
    motif_domain = {}
    for line in motif_domain_table:
        line = line.strip("\n").split("\t")
        print(line)
        if line[0] not in motif_domain:
            motif_domain[line[0]] = []
        motif_domain[line[0]].append(line[1])


#pfam(key) - uniprot(value list) dictionary
with open ("/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/BT_OMV/resources/OMV/extended_OMV_list/OMV_extended_domains.tsv", "r") as protein_domain:
    pfam_uniprot = {}
    for line in protein_domain:
        line = line.strip().split("\t")
        if line[1] not in pfam_uniprot:
            pfam_uniprot[line[1]] = []
        pfam_uniprot[line[1]].append(line[0])


#uniprot(key) - motif(value list) dictionary
uniprot_motif = {}
for key in human.keys():
    for motif in elm_regex:
        match = re.search(str(elm_regex[motif]), str(human[key]))
        if match:
            if rename(key) not in uniprot_motif:
                uniprot_motif[rename(key)] = []
            uniprot_motif[rename(key)].append((motif,str(match.start()),str(match.end())))


with open ("/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/BT_OMV/result/healthy_data/HMI_OMV_extended/macrophage_HMI_H.tsv", "w") as output:
    for pfam, uniprot_list in pfam_uniprot.items():
        for uniprot in uniprot_list:
            for motif in motif_domain:
                if pfam in motif_domain[motif]:
                    for uni, motif_list in uniprot_motif.items():
                        for motif_2 in motif_list:
                            if motif_2[0] == motif:
                                output.write(uni + ";" + ";".join(motif_2) + ";" + ";" + pfam+ ";" + uniprot + "\n")
