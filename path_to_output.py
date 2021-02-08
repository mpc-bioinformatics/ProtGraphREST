
from sanic import exceptions, response
import igraph
import os

from typing import List, Union, Optional, Literal
import falcon 
import json

from pydantic import BaseModel, ValidationError, validator
def parse_path(path: str):


    return None




class Path(BaseModel):
    # Parse it first to List of ints, then strings (comma or -> seperated)
    path: Optional[List[int]] = []
    # Parse a list of strings (comma or -> seperated), then strings (; sperated (comma or -> seperated))
    paths: Optional[List[List[int]]] = []
    returns: Optional[Literal["text", "json"]] = "text"

    @validator("path", pre=True)
    def convert_str_to_list(cls, v):
        if type(v) == str:
            return v.replace("->", ",").split(",")
        return v

    @validator("paths", pre=True)
    def convert_str_to_list_of_lists(cls, v):
        if type(v) == str:
            return [x.replace("->", ",").split(",") for x in v.split(";")]
        elif type(v) == list:
            for idx, x in enumerate(v):
                if type(x) == str:
                    v[idx] = x.replace("->", ",").split(",")
        return v




class PathToPeptide(object):

    def __init__(self, base_dir):
        self.base_dir = base_dir

    def on_get(self, req, resp, accession):
        try:
            # Parse only the Query Parameters
            path_obj = Path(**req.params)
        except ValidationError as ve:
            # TODO raise proper exception!
            resp.body = ve.json()
            resp.status = falcon.HTTP_400
            return

        paths = [] 
        paths.extend(path_obj.paths)
        if len(path_obj.path) != 0:
            paths.append(path_obj.path)
        if len(paths) == 0:
            # TODO throw proper exception
            resp.body = "A path needs to be provided"
            resp.status = falcon.HTTP_400
            return 

        # Check if alphanumeric
        if not accession.isalnum():
            # TODO throw proper exception
            resp.body = "Accession can only consist of +[a-zA-Z0-9]"
            resp.status = falcon.HTTP_400
            return 

        # Get directory (non flat structure) # TODO maybe we should allow both: flat and nonflat?
        prot_graph_path = os.path.join(
            self.base_dir,
            *[x for x in accession[:-1]],
            accession[-1] + ".pickle"
        )

        if not os.path.isfile(prot_graph_path):
            # TODO throw proper exception
            resp.body = "Graph of Protein does not exist!"
            resp.status = falcon.HTTP_404
            return 
        # check if file exists



        # Load graph
        graph = igraph.read(prot_graph_path)

        # For each path retrieve the peptide sequence:
        peptides = []
        for path in paths:
            if max(path) > graph.vcount():
                # TODO throw proper exception
                resp.body = "Path {} has invalid vertex ids (too large)".format(path)
                resp.status = falcon.HTTP_400
                return 

            # Check if the path is connected
            if not all(map(lambda x: graph.are_connected(x[0], x[1]), zip(path, path[1:]))):
                # TODO throw proper exception
                resp.body = "Path {} is not connected".format(path)
                resp.status = falcon.HTTP_400
                return 

            # Check if path goes from start to end:
            if graph.vs[path[0]]["aminoacid"] != "__start__" or \
            graph.vs[path[-1]]["aminoacid"] != "__end__":
                # TODO throw proper exception
                resp.body = "Path {} does not go from the start node to the end node".format(path)
                resp.status = falcon.HTTP_400
                return

            # Append the peptide to the list
            peptides.append(
                # Strip specific start and end
                "".join(graph.vs[path[1:-1]]["aminoacid"])
            )


        # Return the content depending on return type
        if path_obj.returns == "json":
            resp.set_header("content-type", "application/json")
            resp.body = json.dumps(peptides, ensure_ascii=False)
        else:
            resp.set_header("content-type", "text/plain")
            resp.body = "\n".join(peptides)
        resp.status = falcon.HTTP_200





    def on_post(self, req, resp, accession):
        if int(req.headers["CONTENT-LENGTH"]) != 0 and "CONTENT-TYPE" not in req.headers:
            # TODO throw proper exception
            resp.body = "Content-Type needs to be set to application/json"
            resp.status = falcon.HTTP_400
            return

        if  req.headers["CONTENT-TYPE"] != "application/json":
            # TODO throw proper exception
            resp.body = "Content-Type needs to be set to application/json"
            resp.status = falcon.HTTP_400
            return

        paras = req.media
        try:
            # Parse only the Query Parameters
            path_obj_query = Path(**req.params)
            path_obj_body = Path(**paras)
        except ValidationError as ve:
            # TODO raise proper exception!
            resp.body = ve.json()
            resp.status = falcon.HTTP_400
            return

        paths = [] 
        paths.extend(path_obj_query.paths)
        if len(path_obj_query.path) != 0:
            paths.append(path_obj_query.path)

        paths.extend(path_obj_body.paths)
        if len(path_obj_body.path) != 0:
            paths.append(path_obj_body.path) 
        if len(paths) == 0:
            # TODO throw proper exception
            resp.body = "A path needs to be provided"
            resp.status = falcon.HTTP_400
            return 

        # Check if alphanumeric
        if not accession.isalnum():
            # TODO throw proper exception
            resp.body = "Accession can only consist of +[a-zA-Z0-9]"
            resp.status = falcon.HTTP_400
            return 

        # Get directory (non flat structure) # TODO maybe we should allow both: flat and nonflat?
        prot_graph_path = os.path.join(
            self.base_dir,
            *[x for x in accession[:-1]],
            accession[-1] + ".pickle"
        )

        if not os.path.isfile(prot_graph_path):
            # TODO throw proper exception
            resp.body = "Graph of Protein does not exist!"
            resp.status = falcon.HTTP_404
            return 
        # check if file exists

        # Load graph
        graph = igraph.read(prot_graph_path)

        # For each path retrieve the peptide sequence:
        peptides = []
        for path in paths:
            if max(path) > graph.vcount():
                # TODO throw proper exception
                resp.body = "Path {} has invalid vertex ids (too large)".format(path)
                resp.status = falcon.HTTP_400
                return 

            # Check if the path is connected
            if not all(map(lambda x: graph.are_connected(x[0], x[1]), zip(path, path[1:]))):
                # TODO throw proper exception
                resp.body = "Path {} is not connected".format(path)
                resp.status = falcon.HTTP_400
                return 

            # Check if path goes from start to end:
            if graph.vs[path[0]]["aminoacid"] != "__start__" or \
            graph.vs[path[-1]]["aminoacid"] != "__end__":
                # TODO throw proper exception
                resp.body = "Path {} does not go from the start node to the end node".format(path)
                resp.status = falcon.HTTP_400
                return

            # Append the peptide to the list
            peptides.append(
                # Strip specific start and end
                "".join(graph.vs[path[1:-1]]["aminoacid"])
            )


        # Return the content depending on return type
        if "json" in [path_obj_query.returns, path_obj_body.returns]:
            resp.set_header("content-type", "application/json")
            resp.body = json.dumps(peptides, ensure_ascii=False)
        else:
            resp.set_header("content-type", "text/plain")
            resp.body = "\n".join(peptides)
        resp.status = falcon.HTTP_200




class PathToFasta(object):

    def __init__(self, base_dir):
        self.base_dir = base_dir

    def on_get(self, req, resp, accession):
        try:
            # Parse only the Query Parameters
            path_obj = Path(**req.params)
        except ValidationError as ve:
            # TODO raise proper exception!
            resp.body = ve.json()
            resp.status = falcon.HTTP_400
            return

        paths = [] 
        paths.extend(path_obj.paths)
        if len(path_obj.path) != 0:
            paths.append(path_obj.path)
        if len(paths) == 0:
            # TODO throw proper exception
            resp.body = "A path needs to be provided"
            resp.status = falcon.HTTP_400
            return 

        # Check if alphanumeric
        if not accession.isalnum():
            # TODO throw proper exception
            resp.body = "Accession can only consist of +[a-zA-Z0-9]"
            resp.status = falcon.HTTP_400
            return 

        # Get directory (non flat structure) # TODO maybe we should allow both: flat and nonflat?
        prot_graph_path = os.path.join(
            self.base_dir,
            *[x for x in accession[:-1]],
            accession[-1] + ".pickle"
        )

        if not os.path.isfile(prot_graph_path):
            # TODO throw proper exception
            resp.body = "Graph of Protein does not exist!"
            resp.status = falcon.HTTP_404
            return 
        # check if file exists



        # Load graph
        graph = igraph.read(prot_graph_path)

        # For each path retrieve the peptide sequence:
        peptides = []
        for path in paths:
            if max(path) > graph.vcount():
                # TODO throw proper exception
                resp.body = "Path {} has invalid vertex ids (too large)".format(path)
                resp.status = falcon.HTTP_400
                return 

            # Check if the path is connected
            if not all(map(lambda x: graph.are_connected(x[0], x[1]), zip(path, path[1:]))):
                # TODO throw proper exception
                resp.body = "Path {} is not connected".format(path)
                resp.status = falcon.HTTP_400
                return 

            qualifiers = []
            for x, y in zip(path, path[1:]):
                edge = graph.es.find(_between=((x,), (y,)))

                if "qualifiers" in edge.attributes():
                    for q in edge["qualifiers"]:
                        qualifiers.append(q.type)

            # Check if path goes from start to end:
            if graph.vs[path[0]]["aminoacid"] != "__start__" or \
            graph.vs[path[-1]]["aminoacid"] != "__end__":
                # TODO throw proper exception
                resp.body = "Path {} does not go from the start node to the end node".format(path)
                resp.status = falcon.HTTP_400
                return

            # Append the peptide to the list
            peptides.append(
                # Strip specific start and end
               (
                   "".join(graph.vs[path[1:-1]]["aminoacid"]),
                    qualifiers
               )
            )


        # Return the content depending on return type
        if path_obj.returns == "json":
            content = []
            for p in peptides:
                content.append( 
                    {
                        "head" : ">lcl|PEPTIDE_" + accession + "|PATH=" + "->".join([str(i) for i in path]) + "|QUALIFIERS=" + ",".join(p[1]),
                        "seq" : p[0]
                    }
                )
            resp.set_header("content-type", "application/json")
            resp.body = json.dumps(content, ensure_ascii=False)
        else:
            content = ""
            for p in peptides:
                content += ">lcl|PEPTIDE_" + accession + "|PATH=" + "->".join([str(i) for i in path]) + "|QUALIFIERS=" + ",".join(p[1])
                content += "\n" +  '\n'.join(p[0][i:i+3] for i in range(0, len(p[0]), 60))  + "\n"
            resp.set_header("content-type", "text/plain")
            resp.body = content
        resp.status = falcon.HTTP_200





    def on_post(self, req, resp, accession):
        if int(req.headers["CONTENT-LENGTH"]) != 0 and "CONTENT-TYPE" not in req.headers:
            # TODO throw proper exception
            resp.body = "Content-Type needs to be set to application/json"
            resp.status = falcon.HTTP_400
            return

        if  req.headers["CONTENT-TYPE"] != "application/json":
            # TODO throw proper exception
            resp.body = "Content-Type needs to be set to application/json"
            resp.status = falcon.HTTP_400
            return

        paras = req.media
        try:
            # Parse only the Query Parameters
            path_obj_query = Path(**req.params)
            path_obj_body = Path(**paras)
        except ValidationError as ve:
            # TODO raise proper exception!
            resp.body = ve.json()
            resp.status = falcon.HTTP_400
            return

        paths = [] 
        paths.extend(path_obj_query.paths)
        if len(path_obj_query.path) != 0:
            paths.append(path_obj_query.path)

        paths.extend(path_obj_body.paths)
        if len(path_obj_body.path) != 0:
            paths.append(path_obj_body.path) 
        if len(paths) == 0:
            # TODO throw proper exception
            resp.body = "A path needs to be provided"
            resp.status = falcon.HTTP_400
            return 

        # Check if alphanumeric
        if not accession.isalnum():
            # TODO throw proper exception
            resp.body = "Accession can only consist of +[a-zA-Z0-9]"
            resp.status = falcon.HTTP_400
            return 

        # Get directory (non flat structure) # TODO maybe we should allow both: flat and nonflat?
        prot_graph_path = os.path.join(
            self.base_dir,
            *[x for x in accession[:-1]],
            accession[-1] + ".pickle"
        )

        if not os.path.isfile(prot_graph_path):
            # TODO throw proper exception
            resp.body = "Graph of Protein does not exist!"
            resp.status = falcon.HTTP_404
            return 
        # check if file exists

        # Load graph
        graph = igraph.read(prot_graph_path)

        # For each path retrieve the peptide sequence:
        peptides = []
        for path in paths:
            if max(path) > graph.vcount():
                # TODO throw proper exception
                resp.body = "Path {} has invalid vertex ids (too large)".format(path)
                resp.status = falcon.HTTP_400
                return 

            # Check if the path is connected
            if not all(map(lambda x: graph.are_connected(x[0], x[1]), zip(path, path[1:]))):
                # TODO throw proper exception
                resp.body = "Path {} is not connected".format(path)
                resp.status = falcon.HTTP_400
                return 

            qualifiers = []
            for x, y in zip(path, path[1:]):
                edge = graph.es.find(_between=((x,), (y,)))

                if "qualifiers" in edge.attributes():
                    for q in edge["qualifiers"]:
                        qualifiers.append(q.type)

            # Check if path goes from start to end:
            if graph.vs[path[0]]["aminoacid"] != "__start__" or \
            graph.vs[path[-1]]["aminoacid"] != "__end__":
                # TODO throw proper exception
                resp.body = "Path {} does not go from the start node to the end node".format(path)
                resp.status = falcon.HTTP_400
                return

            # Append the peptide to the list
            peptides.append(
                # Strip specific start and end
               (
                   "".join(graph.vs[path[1:-1]]["aminoacid"]),
                    qualifiers
               )
            )


        # Return the content depending on return type
        if "json" in [path_obj_query.returns, path_obj_body.returns]:
            content = []
            for p in peptides:
                content.append( 
                    {
                        "head" : ">lcl|PEPTIDE_" + accession + "|PATH=" + "->".join([str(i) for i in path]) + "|QUALIFIERS=" + ",".join(p[1]),
                        "seq" : p[0]
                    }
                )
            resp.set_header("content-type", "application/json")
            resp.body = json.dumps(content, ensure_ascii=False)
        else:
            content = ""
            for p in peptides:
                content += ">lcl|PEPTIDE_" + accession + "|PATH=" + "->".join([str(i) for i in path]) + "|QUALIFIERS=" + ",".join(p[1])
                content += "\n" +  '\n'.join(p[0][i:i+3] for i in range(0, len(p[0]), 60))  + "\n"
            resp.set_header("content-type", "text/plain")
            resp.body = content
        resp.status = falcon.HTTP_200

