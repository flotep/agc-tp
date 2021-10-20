#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
import operator

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as file:
        sequence = ""
        for line in file :
            if not line.startswith(">"):
                sequence = sequence + line.strip()  
            else:
                if len(sequence) >= minseqlen :
                    yield sequence
                sequence = ""
        yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    sequences = read_fasta(amplicon_file, minseqlen)
    counter = Counter(sequences)
    dict = {k: v for k, v in counter.items() if v >=  mincount}
    dict_sorted = sorted(dict.items(), key=operator.itemgetter(1), reverse=True)
    
    for seq,count in dict_sorted :
        yield [seq,count]

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size) :
    seq_list = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    kmer_dict = {}
    kmer_dict = get_unique_kmer(kmer_dict, seq_list[0][0], 0, kmer_size)
    yield seq_list[0]
    kmer_dict = get_unique_kmer(kmer_dict, seq_list[1][0], 1, kmer_size)
    yield seq_list[1]

    for i, sequence in enumerate (seq_list):
        if i > 1 : 
            perc_identity_matrix = []
            chunks = get_chunks(sequence[0], chunk_size)
            for chunk in range (len(chunks)) : 
                
                best_mates = search_mates(kmer_dict, chunks[chunk], kmer_size)
            
                parent1 = get_chunks(seq_list[best_mates[0]][0],chunk_size)
                parent2 = get_chunks(seq_list[best_mates[1]][0],chunk_size)

                identity_1 = get_identity(nw.global_align(chunks[chunk], parent1[chunk], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH"))))
                identity_2 = get_identity(nw.global_align(chunks[chunk], parent2[chunk], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH"))))
                perc_identity_matrix.append([identity_1, identity_2])
            if detect_chimera(perc_identity_matrix) == False :
                kmer_dict = get_unique_kmer(kmer_dict, sequence[0], i, kmer_size)
                yield([sequence[0], sequence[1]])

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    otu_list = []
    seq_list = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    

    #Pour chaque sequence [i] dans la liste de séquence
    for i, sequence  in enumerate(seq_list):
        if i == 0 : #la 1ère séquence est une OTU
            otu_list.append(sequence)

        else : #on aligne la séquence [i] avec les OTUs [j]
            identity = False
            for j in range (len(otu_list)):
                align = nw.global_align(seq_list[i][0], otu_list[j][0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
                
                if get_identity(align) > 97 and otu_list[j][1] > seq_list[i][1] : #si on ne trouve pas de séquence identique et plus abondante dans la liste d'OTUs
                    identity = True
                    break

            if identity == False :
                otu_list.append(sequence) #cette séquence est un OTU 

    return otu_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as file:
        for i, OTU in enumerate(OTU_list):
            file.write(">OTU_{} occurrence:{}\n".format(i+1, OTU[1]))
            file.write(fill(OTU[0]) + "\n")

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size): #OK
    kmer_list = list(cut_kmer(sequence, kmer_size))
    for kmer in kmer_list:
        if kmer in kmer_dict: #si le kmer est dejà dans le dict
            kmer_dict[kmer] += [id_seq] #on ajoute l'id de la séquence aux autres id 
        else:
            kmer_dict[kmer] = [id_seq] #sinon on initialise ce nouveau kmer dans le dico
    return kmer_dict
    
def search_mates(kmer_dict, sequence, kmer_size):
    best_mates = []
    id_seqs = []
    kmer_list = list(cut_kmer(sequence, kmer_size)) #on coupe le chunk en kmers
    
    for kmer in kmer_list:
        if kmer in kmer_dict.keys(): #si un kmer du chunk correspond à un kmer recensé dans le dict
            id_seqs = id_seqs + kmer_dict[kmer] #on récupère les id des séquences correspondant à ce kmer dans le dict
    
    #Les 2 sequences id qui présentent le plus de kmer en commun avec le chunk sont les parents
    parents = Counter(id_seqs).most_common(2) 
    for parent in parents : 
        best_mates.append(parent[0])

    return best_mates

def std(data):
    return statistics.stdev(data)

def mean(data):
    return statistics.mean(data)

def detect_chimera(perc_identity_matrix): #OK
    similarity_seq1 = False
    similarity_seq2 = False

    #Conditions pour retourner True :
    #A) Moyenne des écarts types > 5 
    if mean([std(perc) for perc in perc_identity_matrix]) <= 5 :
        return False
    else : 
        while similarity_seq1 == False and similarity_seq2 == False :
            for perc in perc_identity_matrix : 
                if perc[0] > perc[1]:
                    similarity_seq1 = True
                elif perc[1] > perc[0]:
                    similarity_seq2 = True

        #B) Au moins 1 segment plus similaire à la seq1 et un autre plus similaire à la seq2
        if similarity_seq1 == True and similarity_seq2 == True :
            return True 
        else : 
            return False


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()
