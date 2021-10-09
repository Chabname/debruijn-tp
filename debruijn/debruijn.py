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

"""Perform assembly based on debruijn graph."""

import argparse
import os
#import sys
#import matplotlib
#from operator import itemgetter
import random
from random import randint
#from numpy import true_divide
import statistics
import networkx as nx
random.seed(9001)


__author__ = "Chabname GHASSEMI"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Chabname GHASSEMI"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Chabname GHASSEMI"
__email__ = "chabname@wanadoo.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{path} is a directory"
        else:
            msg = "{path} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{sys.argv[0]} -h")
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """
    Reading fastq file and create sequence generator
    Parameters:
        fastq_file: fasta file to read
    Returns:
        sequence (yield is the output)
    """
    path = isfile(fastq_file)
    with open(path, "r", encoding="utf-8") as file:
        for _ in file:
            try:
                yield next(file).strip()
                next(file)
                next(file)
            except StopIteration:
                return



def cut_kmer(read, kmer_size):
    """
    Getting kmer for a sequence
    Parameters:
        read: a sequence (read) we want to cut
        kmer_size: size of the kmer (length of cut)
    Returns:
        dictionnary: a dictionnary of kmer with their occurences
    """
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """
    Building the dictionnary with number of occurences
    Parameters:
        fastq_file: fasta file
        kmer_size: size of the kmer
    Returns:
        dictionnary: a dictionnary of kmer with their occurences
    """
    dictionnary = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in dictionnary.keys():
                dictionnary[kmer] += 1
            else:
                dictionnary[kmer] = 1
    return dictionnary



def build_graph(kmer_dict):
    """
    Building a weighted directed edge graph
    Parameters:
        kmer_dict: dictionnary of kmer with their occurences
    Returns:
        diGraph: a weighted directed edge graph
    """
    digraph = nx.DiGraph()

    for key_dic, dic_weight in kmer_dict.items():
        # key_dic[:-1] : get the kmer without the last amino acid
        # key_dic[1:] : get the kmer without the first amino acid
        print(key_dic[1:])
        digraph.add_edge(key_dic[:-1], key_dic[1:], weight = dic_weight)
    return digraph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    Remove one or multiple path
    Parameters:
        graph: tree
        path_list: the path to remove
        delete_entry_node: boolean to know if starting node is removed
        delete_sink_node: boolean to know if ending node is removed
    Returns:
        Remove the path
    """
    for path in path_list:
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])

        graph.remove_nodes_from(path[1:-1])
    return graph



def std(data):
    """
    Standard deviation

    Parameters:
        data: List of values
    Returns:
        Return the standard deviation
    """
    return statistics.stdev(data)



def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """
    Select the best and remove the other path
    Parameters:
        graph: the fraph we want to simplify
        path_list: list of path to check
        path_length: lenth of path
        weight_avg_list: list of the average weight of the path (in the same order of the list path)
        delete_entry_node: boolean to know if starting node is removed
        delete_sink_node: boolean to know if ending node is removed
    Returns:
        Return the best path: a graph
    """
    remove_path = path_list
    index_max = 0

    # Best frequency on average
    if std(weight_avg_list) != 0:
        max_value = max(weight_avg_list)
        index_max = weight_avg_list.index(max_value)
    # Best frequency on length
    elif std(path_length) != 0:
        max_value = max(path_length)
        index_max = path_length.index(max_value)
    #random
    else:
        index_max = randint(0, len(path_list))

    remove_path.pop(index_max)

    return remove_paths(graph, remove_path, delete_entry_node, delete_sink_node)




def path_average_weight(graph, path):
    """
    Get the average weight of a path
    Parameters:
        graph: tree
        path: one possible path (contigs)
    Returns:
        The average weight of a path
    """
    list_weight = []
    sub_graphes = graph.subgraph(path).edges(data=True)
    for edge in sub_graphes:
        list_weight.append(edge[2]["weight"])
    return statistics.mean(list_weight)


def solve_bubble(graph, ancestor_node, descendant_node):
    """
    Get the average weight of a path
    Parameters:
        graph: tree
        path: one possible path (contigs)
    Returns:
        The average weight of a path
    """
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = []
    weight_avg_list = []
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))

    return select_best_path(graph, path_list, path_length, weight_avg_list)


def simplify_bubbles(graph):
    """
    Get the average weight of a path
    Parameters:
        graph: tree
        path: one possible path (contigs)
    Returns:
        The average weight of a path
    """
    bubble = False
    node = []
     # For every nodes
    for node in graph.nodes:
        # Check if the node exists (because of the recursivity at the end)
        if node in graph.nodes:
            # get all predecessors
            predecessors = list(graph.predecessors(node))
            if len(predecessors) > 1:
                for i, pred in enumerate(predecessors):
                    for j in range(i+1, len(predecessors)):
                        # Get the nearest ancestor between predecessor i eand predessor j
                        ancestor = nx.lowest_common_ancestor(graph, pred, predecessors[j])
                        if ancestor is not None:
                            bubble = True
                            break

    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor, node))

    return graph



def solve_entry_tips(graph, starting_nodes):
    """
    Solve the multiple start nodes
    Parameters:
        graph: tree
        starting_nodes: List of all possible starting nodes
    Returns:
        Cleaned graph to have only one "entry" path
    """
    path_list = []
    # We have to find all nodes having mulitiple ancestor
    # For every nodes
    for node in graph.nodes:
        # Check if the node exists
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for start in starting_nodes:
                #This is not a bubble so there is only one path
                path_list += list(nx.all_simple_paths(graph, start, node))
    #Then we check if we have nodes with multiple predessors
    if len(path_list) != 0:
        path_length = []
        weight_avg_list = []
        for path in path_list:
            path_length.append(len(path))
            weight_avg_list.append(path_average_weight(graph, path))

    #we have to remove the first node and all betWeen but not the last node
    return select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=True)


def solve_out_tips(graph, ending_nodes):
    """
    Solve the multiple end nodes
    Parameters:
        graph: tree
        starting_nodes: List of all possible starting nodes
    Returns:
        Cleaned graph to have only one "entry" path
    """
    path_list = []
    # We have to find all nodes having mulitiple ancestor
    # For every nodes
    for node in graph.nodes:
        # Check if the node exists
        successors = list(graph.successors(node))
        if len(successors) > 1:
            for end in ending_nodes:
                #This is not a bubble so there is only one path
                path_list += list(nx.all_simple_paths(graph, node, end))

    #Then we check if we have nodes with multiple predessors
    if len(path_list) != 0:
        path_length = []
        weight_avg_list = []
        for path in path_list:
            path_length.append(len(path))
            weight_avg_list.append(path_average_weight(graph, path))

    #we have to remove the last node and all betWeen but not the first node
    return select_best_path(graph, path_list, path_length, weight_avg_list, delete_sink_node=True)


def get_starting_nodes(graph):
    """
    Find all nodes which don't have predessors
    Parameters:
        graph: tree of prefixes and suffixes kmers
    Returns:
        A list of starting nodes
    """
    start_node = []
    for node in graph:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) == 0:
            start_node.append(node)
    return start_node



def get_sink_nodes(graph):
    """
    Fin all nodes which don't have successors
    Parameters:
        graph: tree of prefixes and suffixes kmers
    Returns:
        A list of final nodes
    """
    sink_node = []
    for node in graph:
        successsors = list(graph.successors(node))
        if len(successsors) == 0:
            sink_node.append(node)
    return sink_node



def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Get all the possible contiges with simple path
    Parameters:
        graph:tree of prefixes and suffixes kmers
    	starting_nodes: list of starting nodes
    	ending_nodes: list of end nodes
    Returns:
        A list of tuples with contig and length of contig
    """
    contigs=[]
    for start_node in starting_nodes:
        for end_node in ending_nodes:

            if nx.has_path(graph, start_node, end_node):
                constructed_contig = ""
                for path in nx.all_simple_paths(graph, start_node, end_node):

                    for index, word in enumerate(path):
                        if index == 0 :
                            constructed_contig = word
                        else:
                            constructed_contig += word[-1]
                    contigs.append([constructed_contig, len(constructed_contig)])
    return contigs



def save_contigs(contigs_list, output_file):
    """
    Save the finded contigs and write them into a file
    Parameters:
        contigs_list: List of all contigs with their length
    	output_file: path of the output file
    Returns:
        A list of tuples with contig and length of contig
    """
    with open(output_file, "w", encoding="utf-8") as file:
        for index, contig in enumerate(contigs_list):
            file.write(">contig_" + str(index) + " len=" + str(contig[1]) + "\n")
            file.write(fill(contig[0]) + "\n")



def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    #fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    #plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    #with open(graph_file, "wt") as save:
    #    pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_dictionnary = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dictionnary)
    graph = simplify_bubbles(graph)
    sart_nodes = get_starting_nodes(graph)
    end_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, sart_nodes)
    graph = solve_out_tips(graph, end_nodes)
    contigs = get_contigs(graph, sart_nodes, end_nodes)
    save_contigs(contigs, "data/results.txt")

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
