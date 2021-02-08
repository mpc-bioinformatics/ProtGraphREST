import argparse
import os

from typing import List
import igraph
# from sanic import Sanic, exceptions, response
# from sanic_openapi import doc, swagger_blueprint

import path_to_output
import get_parameters as get_p
import models
from fastapi.exceptions import RequestValidationError


from waitress import serve
import json

import falcon
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


if __name__ == '__main__':
    GLOABL_ARGS = parse_args()

    app.add_route("/{accession}/path_to_peptide", path_to_output.PathToPeptide(GLOABL_ARGS["base_folder"]))
    app.add_route("/{accession}/path_to_fasta", path_to_output.PathToFasta(GLOABL_ARGS["base_folder"]))

    serve(app, listen="*:8000")

    # Start Sanic
    # uvicorn.run(app, host="0.0.0.0", port=8000)
