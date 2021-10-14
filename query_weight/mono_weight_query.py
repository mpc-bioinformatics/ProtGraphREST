import ctypes
import inspect
import json
import threading
import time
import timeit

import falcon
import igraph
import numpy as np
from protgraph.aa_masses_annotation import _get_mass_dict

import query_weight.query_algorithms as qa
from graph_utils import get_graph_path, get_pdb_path, get_start_and_end_node
from models import MonoWeigthQuery
from models_utils import load_model
from prot_graph_exception import ProtGraphException

ALGORITHMS = dict(
        top_sort=qa.top_sort_query,
        bfs_fifo=qa.bfs_fifo,
        bfs_filo=qa.bfs_fifo,
        dfs=qa.dfs,
        top_sort_attrs=qa.top_sort_attrs_query,
        top_sort_attrs_limit_var=qa.top_sort_attrs_query_limit_variants
    )


def _check_header(req):
    """ checks if the header for POST is set """
    if int(req.headers["CONTENT-LENGTH"]) != 0 and "CONTENT-TYPE" not in req.headers:
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Content-Length needs to be set to application/json"}, indent=4)
        )
    if req.headers["CONTENT-TYPE"] != "application/json":
        raise ProtGraphException(
            falcon.HTTP_400,
            json.dumps({"message": "Content-Type needs to be set to application/json"}, indent=4)
        )


class QueryWeight(object):
    """ This is a generic class, utilizing the method, which generate paths from qeight queries. """

    def __init__(self, base_dir, weight_factor, method):
        self.base_dir = base_dir
        self.weight_factor = weight_factor
        self.method = method
        self.mass_dict = _get_mass_dict(factor=self.weight_factor)

    def _return_content(self, resp, peptides, peptide_weights, peptide_seqs, time):
        """ Return the content. Speicifically Peptide -Path, -Weight and -Sequence. """
        # Generate returning dict
        return_dict = dict(
            time=time
        )
        # Returning the paths for the peptide
        return_dict["results"] = [
            dict(path=p, weight=w, seq=s)
            for p, w, s in zip(peptides, peptide_weights, peptide_seqs)
        ]

        resp.set_header("content-type", "application/json")
        resp.body = json.dumps(return_dict, ensure_ascii=False)

    def _execute_query(self, query, accession):
        """ executes the weight query and gathers other information """
        # Get Protein depending on accession and path
        prot_graph_path = get_graph_path(self.base_dir, accession)

        # Load graph
        graph = igraph.read(prot_graph_path)

        # Get pdb if not generate it and get it!
        n_pdb = get_pdb_path(self.base_dir, accession, graph, query.k)

        # Get intervals
        w = self.weight_factor*query.mono_weight
        if query.unit == "Da":
            wf = self.weight_factor*query.mass_tolerance
            q_interval = np.array([
                (w - wf),
                (w + wf),
            ])
        else:  # == "ppm"
            q_interval = np.array([
                w - (w / 1000000) * query.mass_tolerance,
                w + (w / 1000000) * query.mass_tolerance
            ])

        # Set start and end node
        start, end = get_start_and_end_node(graph)

        # Execute and measure time  TODO measure time
        resulting_paths = []

        def _async_raise(tid, exctype):
            '''Raises an exception in the threads with id tid'''
            if not inspect.isclass(exctype):
                raise TypeError("Only types can be raised (not instances)")
            res = ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), ctypes.py_object(exctype))
            if res == 0:
                raise ValueError("invalid thread id")
            elif res != 1:
                # "if it returns a number greater than one, you're in trouble,
                # and you should call it again with exc=NULL to revert the effect"
                ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), None)
                raise SystemError("PyThreadState_SetAsyncExc failed")

        class ThreadWithExc(threading.Thread):

            def __init__(self, method):
                super().__init__()
                self.method = method

            def run(self):
                try:
                    resulting_paths.extend(self.method(start, end, q_interval, graph, n_pdb))
                finally:
                    pass

            def _get_my_tid(self):
                if not self.isAlive():
                    raise threading.ThreadError("the thread is not active")

                # do we have it cached?
                if hasattr(self, "_thread_id"):
                    return self._thread_id

                # no, look for it in the _active dict
                for tid, tobj in threading._active.items():
                    if tobj is self:
                        self._thread_id = tid
                        return tid

                raise AssertionError("could not determine the thread's id")

            def raiseExc(self, exctype):
                _async_raise(self._get_my_tid(), exctype)

        t = ThreadWithExc(self.method)

        # STAR MEASURING HERE!
        starttime = timeit.default_timer()
        t.start()
        t.join(query.timeout)
        time_taken = timeit.default_timer() - starttime
        # STOP MEASURING HERE!

        # Set timeout limit
        if time_taken > query.timeout:
            time_taken = query.timeout

        # Terminate thread if still running
        if t.is_alive():
            while t.isAlive():
                time.sleep(0.01)
                try:
                    t.raiseExc(Exception)
                except Exception:
                    pass
            t.join()

        # Get weights, which were actually retrieved:  (via protgraph)
        resulting_weights = [
            sum(
                [
                    self.mass_dict[x][0]
                    for x in "".join(graph.vs[path]["aminoacid"])
                    .replace("__start__", "")
                    .replace("__end__", "")
                ]
            ) / self.weight_factor
            for path in resulting_paths
        ]

        # Get weights, which were actually retrieved:  (via protgraph)
        resulting_seq = [
            "".join(graph.vs[path]["aminoacid"])
            .replace("__start__", "")
            .replace("__end__", "")
            for path in resulting_paths
        ]

        return resulting_paths, resulting_weights, resulting_seq, time_taken

    def on_get(self, req, resp, accession):
        # Load Query
        query = load_model(MonoWeigthQuery, req.params)

        # Get peptides
        results = self._execute_query(query, accession)

        # Return the content depending on return type
        self._return_content(resp, *results)

    def on_post(self, req, resp, accession):
        # Check headers
        _check_header(req)

        # Load Query
        query = load_model(MonoWeigthQuery, req.media)

        # Get peptides
        results = self._execute_query(query, accession)

        # Return the content depending on return type
        self._return_content(resp, *results)
