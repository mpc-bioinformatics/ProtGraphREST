import json
import os

import falcon

from prot_graph_exception import ProtGraphException


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
