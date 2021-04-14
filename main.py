import argparse
import os

import falcon
from falcon_swagger_ui import register_swaggerui_app
from waitress import serve

import path_to_output
from prot_graph_exception import ProtGraphException
from query_weight import mono_weight_query as wq

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
    parser.add_argument(
        "--mass_dict_factor", "-mdf", type=float, default=1000000000,
        help="Set the factor for the masses which was used to generate the graphs. "
        "The default is set to 1 000 000 000, so that each mass can be converted into integers."
    )

    args = parser.parse_args()

    return dict(
        base_folder=args.base_folder,
        mass_dict_factor=args.mass_dict_factor
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

    # Add resources folder to:
    app.add_static_route('/resources', os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "resources"
        )
    )

    register_swaggerui_app(
        app, "/swagger", "/resources/openapi.yaml",
        page_title="Swagger for ProtGraphREST"
    )

    app.add_route("/{accession}/path_to_peptide", path_to_output.PathToPeptide(GLOABL_ARGS["base_folder"]))
    app.add_route("/{accession}/path_to_fasta", path_to_output.PathToFasta(GLOABL_ARGS["base_folder"]))

    # Routes for weight queries
    app.add_route(
        "/{accession}/top_sort/query_mono_weight",
        wq.QueryWeight(GLOABL_ARGS["base_folder"], GLOABL_ARGS["mass_dict_factor"], wq.ALGORITHMS["top_sort"])
    )
    app.add_route(
        "/{accession}/bfs_fifo/query_mono_weight",
        wq.QueryWeight(GLOABL_ARGS["base_folder"], GLOABL_ARGS["mass_dict_factor"], wq.ALGORITHMS["bfs_fifo"])
    )
    app.add_route(
        "/{accession}/bfs_filo/query_mono_weight",
        wq.QueryWeight(GLOABL_ARGS["base_folder"], GLOABL_ARGS["mass_dict_factor"], wq.ALGORITHMS["bfs_filo"])
    )
    app.add_route(
        "/{accession}/dfs/query_mono_weight",
        wq.QueryWeight(GLOABL_ARGS["base_folder"], GLOABL_ARGS["mass_dict_factor"], wq.ALGORITHMS["dfs"])
    )

    # Example call for a query via weight
    # http://localhost:8000/A0A4S5AXF8/top_sort/query_mono_weight?unit=ppm&mono_weight=3394.719&mass_tolerance=5
    # http://localhost:8000/P04637/top_sort/query_mono_weight?unit=ppm&mono_weight=1000.719&mass_tolerance=5&timeout=10

    # Example call for getting a peptide:
    # http://localhost:8000/A0A4S5AXF8/path_to_fasta?path=0,24,25,9
    serve(app, listen="*:8000")
