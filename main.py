import argparse
import os

import igraph
from sanic import Sanic, exceptions, response

app = Sanic("ProtGraphREST")


@app.route("<accession:string>/path_to_pep")
async def path_to_peptide(request, accession):
    """
    We load a pickle file and return the peptide from the corresponding protein

    Route Arguments:
    accession -> The UniProt-Accesion of the protein

    ? Arguments:
    path=[List of Integers]

    Returns: The actual peptide as PLAIN text (multiple if more paths are set, seperated by "/n")
    """
    if "path" not in request.args:
        raise exceptions.ServerError("Required Argument: 'path' was not defined!", status_code=400)

    # Get directory (non flat structure) # TODO maybe we should allow both: flat anf nonflat?
    prot_graph_path = os.path.join(
        GLOABL_ARGS["base_folder"],
        *[x for x in accession[:-1]],
        accession[-1] + ".pickle"
    )

    # TODO should we limit accession to [A-Z0-9]?
    if not os.path.isfile(prot_graph_path):
        raise exceptions.ServerError("No graph found for accession {}".format(accession), status_code=404)
    # check if file exists

    # Load graph
    graph = igraph.read(prot_graph_path)

    # For each path retrieve the peptide sequence:
    peptides = []
    for path in request.args["path"]:
        # Parse the integer list
        try:
            path_ints = list(map(int, path.split(",")))
        except Exception:
            raise exceptions.ServerError("Path {} can only consist of ',' and [0-9]", status_code=400)

        # Check if the path is connected
        if not all(map(lambda x: graph.are_connected(x[0], x[1]), zip(path_ints, path_ints[1:]))):
            raise exceptions.ServerError("Path {} is not connected".format(path), status_code=400)

        # Check if path goes from start to end:
        if graph.vs[path_ints[0]]["aminoacid"] != "__start__" or \
           graph.vs[path_ints[-1]]["aminoacid"] != "__end__":
            raise exceptions.ServerError(
                "Path {} does not go from the start node to the end node".format(path), status_code=400
            )

        # Append the peptide to the list
        peptides.append(
            # Strip specific start and end
            "".join(graph.vs[path_ints[1:-1]]["aminoacid"])
        )

    return response.text(
        "\n".join(peptides) + "\n"
    )


def parse_args():
    # Arguments for Parser
    parser = argparse.ArgumentParser(
        description="Graph-Generator-REST-API for generated Graphs of Proteins and Peptides"
    )
    # Needed Arguments for parsing (and additional information for it)
    parser.add_argument(
        "base_folder", type=str,
        help="The base folder containing the generated exported graph files, which can be read by igraph."
    )

    args = parser.parse_args()

    return dict(
        base_folder=args.base_folder
    )


if __name__ == '__main__':
    GLOABL_ARGS = parse_args()

    # Start Sanic
    app.run(host="0.0.0.0", port=8000)
