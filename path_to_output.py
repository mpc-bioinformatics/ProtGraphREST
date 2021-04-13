import json

import falcon
import igraph

from graph_utils import (check_path_incorrect, get_aminoacids, get_graph_path,
                         get_qualifiers)
from models import Path
from models_utils import load_model
from prot_graph_exception import ProtGraphException


def _check_header(req):
    """ checks if the header for POST is set """
    if int(req.headers["CONTENT-LENGTH"]) != 0 and "CONTENT-TYPE" not in req.headers:
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Content-Length needs to be set"}, indent=4)
        )
    if req.headers["CONTENT-TYPE"] != "application/json":
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Content-Type needs to be set to application/json"}, indent=4)
        )


def _concat_paths(*paths: Path):
    """ Concatenate all paths, into one list """
    paths_list = []
    for path in paths:
        paths_list.extend(path.paths)
        if len(path.path) != 0:
            paths_list.append(path.path)
    return paths_list


def _check_paths_length(paths):
    """ Throws an exception if no path is provided. """
    if len(paths) == 0:
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "A path needs to be provided"}, indent=4)
        )


class PathToPeptide(object):

    def __init__(self, base_dir):
        self.base_dir = base_dir

    def _return_content(self, resp, peptides, as_json=True):
        if as_json:
            resp.set_header("content-type", "application/json")
            resp.body = json.dumps(peptides, ensure_ascii=False)
        else:
            resp.set_header("content-type", "text/plain")
            resp.body = "\n".join(peptides)
        resp.status = falcon.HTTP_200

    def _get_peptides(self, resp, prot_graph_path, paths):
        _check_paths_length(paths)

        # Load graph
        graph = igraph.read(prot_graph_path)

        # For each path retrieve the peptide sequence:
        peptides = []
        for path in paths:
            check_path_incorrect(graph, path)
            peptides.append(get_aminoacids(graph, path[1:-1]))

        return peptides

    def on_get(self, req, resp, accession):
        # Get Protein depending on accession and path
        prot_graph_path = get_graph_path(self.base_dir, accession)
        path_obj = load_model(Path, req.params)
        paths = _concat_paths(path_obj)

        # Get peptides
        peptides = self._get_peptides(resp, prot_graph_path, paths)

        # Return the content depending on return type
        self._return_content(resp, peptides, path_obj.returns == "json")

    def on_post(self, req, resp, accession):
        # Check headers
        _check_header(req)
        # Get Protein depending on accession and path
        prot_graph_path = get_graph_path(self.base_dir, accession)
        path_obj_query, path_obj_body = load_model(Path, req.params, req.media)
        paths = _concat_paths(path_obj_query, path_obj_body)

        # Get peptides
        peptides = self._get_peptides(resp, prot_graph_path, paths)

        # Return the content depending on return type
        self._return_content(resp, peptides, "json" in [path_obj_query.returns, path_obj_body.returns])


class PathToFasta(object):

    def __init__(self, base_dir):
        self.base_dir = base_dir

    def _return_content(self, resp, peptides, paths, accession, as_json=True):
        if as_json:
            content = []
            for pep, path in zip(peptides, paths):
                content.append(
                    {
                        "head": ">lcl|PEPTIDE_" + accession
                        + "|PATH=" + "->".join([str(i) for i in path])
                        + "|QUALIFIERS=" + ",".join(pep[1]),
                        "seq": pep[0]
                    }
                )
            resp.set_header("content-type", "application/json")
            resp.body = json.dumps(content, ensure_ascii=False)
        else:
            content = ""
            for pep, path in zip(peptides, paths):
                content += ">lcl|PEPTIDE_" + accession \
                    + "|PATH=" + "->".join([str(i) for i in path]) \
                    + "|QUALIFIERS=" + ",".join(pep[1])
                content += "\n" + '\n'.join(pep[0][i:i+60] for i in range(0, len(pep[0]), 60)) + "\n"
            resp.set_header("content-type", "text/plain")
            resp.body = content
        resp.status = falcon.HTTP_200

    def _get_peptides(self, resp, prot_graph_path, paths):
        _check_paths_length(paths)

        # Load graph
        graph = igraph.read(prot_graph_path)

        # For each path retrieve the peptide sequence:
        peptides = []
        for path in paths:
            check_path_incorrect(graph, path)
            peptides.append(
                (
                    get_aminoacids(graph, path[1:-1]),
                    get_qualifiers(graph, path[1:-1])
                )
            )

        return peptides

    def on_get(self, req, resp, accession):
        # Get Protein depending on accession and path
        prot_graph_path = get_graph_path(self.base_dir, accession)
        path_obj = load_model(Path, req.params)
        paths = _concat_paths(path_obj)

        # Get peptides
        peptides = self._get_peptides(resp, prot_graph_path, paths)

        # Return the content depending on return type
        self._return_content(
            resp, peptides, paths, accession,
            path_obj.returns == "json"
        )

    def on_post(self, req, resp, accession):
        # Check headers
        _check_header(req)
        # Get Protein depending on accession and path
        prot_graph_path = get_graph_path(self.base_dir, accession)
        path_obj_query, path_obj_body = load_model(Path, req.params, req.media)
        paths = _concat_paths(path_obj_query, path_obj_body)

        # Get peptides
        peptides = self._get_peptides(resp, prot_graph_path, paths)

        # Return the content depending on return type
        self._return_content(
            resp, peptides, paths, accession,
            "json" in [path_obj_query.returns, path_obj_body.returns]
        )
