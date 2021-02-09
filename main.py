import argparse

import falcon
from waitress import serve

import path_to_output
from prot_graph_exception import ProtGraphException

app = application = falcon.API()


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


def generic_error_handler(ex, req, resp, params):
    """ Set the generic ErrorHandler for the exception: ProtGraphException"""
    if isinstance(ex, ProtGraphException):
        resp.status = ex.status
        resp.body = ex.body
        for x, y in ex.headers.items():
            resp.set_header(x, y)
    else:
        raise ex


if __name__ == '__main__':
    GLOABL_ARGS = parse_args()

    app.add_error_handler(ProtGraphException, generic_error_handler)

    app.add_route("/{accession}/path_to_peptide", path_to_output.PathToPeptide(GLOABL_ARGS["base_folder"]))
    app.add_route("/{accession}/path_to_fasta", path_to_output.PathToFasta(GLOABL_ARGS["base_folder"]))

    serve(app, listen="*:8000")
