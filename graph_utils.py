import json
import os

import falcon
import numpy as np
from protgraph.export.peptides.pep_fasta import PepFasta

from prot_graph_exception import ProtGraphException
from query_weight.query_algorithms import build_pdb

PF = PepFasta()


def get_pdb_path(base_dir, accession: str, graph, k=5):
    """ Get intervals from graph from numpy file if it exists. If not generate it """
    # Check if accession is correct
    if not accession.isalnum:
        raise ProtGraphException(
            falcon.HTTP_404,
            json.dumps({"message": "Accession can only consist of +[a-zA-Z0-9]"}, indent=4)
        )

    # Get directory (non flat structure) # TODO maybe we should allow both: flat and nonflat?
    path = os.path.join(
        base_dir,
        *[x for x in accession[:-1]],
        accession[-1] + ".pdb" + str(k)
    )

    # Check if the path to the graph file exists
    if not os.path.isfile(path):
        # Then generate pdb
        build_pdb(graph, k=k)

        # Fill entries with nans
        pdb_entries = [x + [[np.nan, np.nan]]*(k - len(x)) for x in graph.vs["pdb"]]

        # Delete entries in graph itself
        del graph.vs["pdb"]

        # Generate numpy matrix out of these intervals
        n_pdb = np.array(pdb_entries)

        # Save it on disk:
        with open(path, "wb") as f:
            np.save(f, n_pdb)

        # Return the interval matrix
        return n_pdb

    # Return the interval matrix
    return np.load(path)


def get_graph_path(base_dir, accession: str):
    """ Gets the path for a protein. Raises exceptions if accession is invalid or if the graph is not existing """
    # Check if accession is correct
    if not accession.isalnum:
        raise ProtGraphException(
            falcon.HTTP_404,
            json.dumps({"message": "Accession can only consist of +[a-zA-Z0-9]"}, indent=4)
        )

    # Get directory (non flat structure) # TODO maybe we should allow both: flat and nonflat?
    path = os.path.join(
        base_dir,
        *[x for x in accession[:-1]],
        accession[-1] + ".pickle"  # TODO only pickle?, do we want to change this?
    )

    # check if the path to the graph file exists
    if not os.path.isfile(path):
        raise ProtGraphException(
            falcon.HTTP_404,
            json.dumps({"message": "Graph of Protein does not exist!"}, indent=4)
        )

    # return if the graph exists
    return path


def check_path_connected(graph, path: list):
    """ raises an exception if path is not connected in graph """
    if not all(map(lambda x: graph.are_connected(x[0], x[1]), zip(path, path[1:]))):
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Path '{}' is not connected".format(path)}, indent=4)
        )


def check_start_end_nodes(graph, path):
    """ check if start and end of a path correspond to the start and end node in a graph """
    if graph.vs[path[0]]["aminoacid"] != "__start__" or \
       graph.vs[path[-1]]["aminoacid"] != "__end__":
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Path {} does not go from the start node to the end node".format(path)}, indent=4)
        )


def check_path_invalid(graph, path):
    """ Checks if a path is invalid, so we check here for connection and if the highest node id is not to large """
    # Check if node_ids are in range
    if max(path) > graph.vcount():
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Path {} has invalid vertex ids (too large).".format(path)}, indent=4)
        )
    if min(path) < 0:
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Path {} cannot have negative vertices.".format(path)}, indent=4)
        )

    # Check if connected
    check_path_connected(graph, path)


def check_path_incorrect(graph, path):
    """
    Checks if a path is incorrect:
    1. Beginning/ending at start/end,
    2. Paths is connected and
    3. The range of number needs to be in range of the graph

    Raises exceptions, if invalid
    """
    check_path_invalid(graph, path)
    check_start_end_nodes(graph, path)


def get_aminoacids(graph, path: list):
    """ Returns the aminoacids attributes (nodes; in order) """

    # return the aminoacids of the path
    return "".join(graph.vs[path]["aminoacid"])


def get_pep_and_header_def(path, graph):
    """ Get the FASTA Header deginition as a string. Also return the peptide """
    peptide = "".join(graph.vs[path[1:-1]]["aminoacid"])
    edges = graph.get_eids(path=path)
    if "cleaved" in graph.es[edges[0]].attributes():
        misses = sum(filter(None, graph.es[edges]["cleaved"]))
    else:
        misses = -1

    acc = PF._get_accession_or_isoform(graph.vs[path[1]])
    start_pos = PF._get_position_or_isoform_position(graph.vs[path[1]])
    end_pos = PF._get_position_or_isoform_position(graph.vs[path[-2]], end=True)
    l_str_qualifiers = PF._get_qualifiers(graph, edges)
    quali_entries = ",".join(l_str_qualifiers)

    part_header = "".join(
        [
            acc, "(", str(start_pos), ":", str(end_pos), ",",
            "mssclvg:", str(misses),
            quali_entries,  ")"
        ]
    )

    return peptide, part_header


def get_start_and_end_node(graph):
    """ Returns the start and end node from a graph (which is unique!) """

    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__end_node__] = graph.vs.select(aminoacid="__end__")
    return __start_node__, __end_node__


def get_qualifiers(graph, path: list):
    """ Returns the qualifiers attributes (edges; in order, if existent) """
    # Get qualifier
    qualifiers = []
    for x, y in zip(path, path[1:]):
        edge = graph.es.find(_between=((x,), (y,)))
        # if it exits
        if "qualifiers" in edge.attributes():
            for q in edge["qualifiers"]:
                qualifiers.append(q.type)

    # return the aminoacids of the path
    return qualifiers
