from collections import defaultdict, deque

import numpy as np

from protgraph.graph_collapse_edges import Or

def _shift_interval_by(intervals, weight):
    """ Shift the intervals by weight """
    return [[x + weight, y + weight] for [x, y] in intervals]


def _merge_overlapping_intervals(intervals):
    """ Get overlapping intervals and merge them """
    intervals = np.array(intervals)
    starts = intervals[:, 0]
    ends = np.maximum.accumulate(intervals[:, 1])
    valid = np.zeros(len(intervals) + 1, dtype=bool)
    valid[0] = True
    valid[-1] = True
    valid[1:-1] = starts[1:] >= ends[:-1]
    return [list(x) for x in np.vstack((starts[:][valid[:-1]], ends[:][valid[1:]])).T]


def _merge_closest_intervals(intervals):
    """ Get the closest interval and merge those two """
    diff = [y[0]-x[1] for x, y in zip(intervals, intervals[1:])]
    argmin = diff.index(min(diff))

    new_interval = [intervals[argmin][0], intervals[argmin+1][1]]
    return intervals[:argmin] + [new_interval] + intervals[argmin+2:]


def build_pdb(graph, k=5):
    """
    Generates the pdb (intervals). Each node will have up to k many intervals
    We build it up via the rev. top. sort. Intervals are merged they overlap.
    If too many intervals are present the closest ones will be merged.

    This uses the mono_weight, but could be extended to use the avrg_weight (TODO DL?)
    """
    rev_top_sort = graph.topological_sorting(mode="IN")

    # Initial attribute values:
    graph.vs[rev_top_sort[0]]["pdb"] = [[0, 0]]

    # iterate
    for node in rev_top_sort[1:]:
        intervals = []
        for out_edge in graph.incident(node, mode="OUT"):
            intervals.extend(
                _shift_interval_by(
                    graph.vs[graph.es[out_edge].target]["pdb"],
                    graph.es[out_edge]["mono_weight"]
                )
            )

        sorted_intervals = _merge_overlapping_intervals(sorted(intervals, key=lambda x: x[0]))

        while True:
            if len(sorted_intervals) <= k:
                break
            else:
                sorted_intervals = _merge_closest_intervals(sorted_intervals)

        graph.vs[node]["pdb"] = sorted_intervals


def dfs(start, stop, tv_interval, _graph, _n_pdb):
    """ Depth First Search of the Graph """
    return _dfs_inner([start.index], stop.index, tv_interval, _graph, _n_pdb)


def _dfs_inner(path, stop, tv, _graph, _n_pdb):
    """ Recursive inner method """
    p = []

    if path[-1] == stop:
        return [path]

    for out_edge in _graph.es.select(_source=path[-1]):
        t_path = path + [out_edge.target]
        edges = _graph.get_eids(path=t_path)
        val_f = _func_dist(_n_pdb[out_edge.target, :], tv - sum(_graph.es[edges]["mono_weight"]))
        if val_f:
            p.extend(_dfs_inner(t_path, stop, tv, _graph, _n_pdb))

    return p


def bfs_filo(start, stop, tv_interval, _graph, _n_pdb):
    """ Breadth-First-Search using a FILO approach. """
    queue = deque()
    queue.append([start.index, [], 0])
    paths = []

    while queue:
        cur_node, how_to_get, sum_weight = queue.pop()  # This is changed!

        if cur_node == stop.index:
            paths.append([*how_to_get, cur_node])
            continue

        edges = _graph.es.select(_source=cur_node)
        eids = [e.index for e in edges]

        targeted_nodes = [e.target for e in edges]
        achieved_tvs = [_graph.es[e]["mono_weight"] + sum_weight for e in eids]
        val_fs = [
            _func_dist(_n_pdb[out_edge.target, :], tv_interval - g_pre)
            for out_edge, g_pre in zip(edges, achieved_tvs)
        ]

        pruned_targeted_nodes = [tn for tn, val_f in zip(targeted_nodes, val_fs) if val_f]
        pruned_achieved_tvs = [tv for tv, val_f in zip(achieved_tvs, val_fs) if val_f]

        for x, y in zip(pruned_targeted_nodes, pruned_achieved_tvs):
            queue.append([x, how_to_get + [cur_node], y])

    return paths


def bfs_fifo(start, stop, tv_interval, _graph, _n_pdb):
    """ Breadth-First-Search using a FIFO approach. (classic approach) """
    queue = deque()
    queue.append([start.index, [], 0])
    paths = []

    while queue:
        cur_node, how_to_get, sum_weight = queue.popleft()

        if cur_node == stop.index:
            paths.append([*how_to_get, cur_node])
            continue

        edges = _graph.es.select(_source=cur_node)
        eids = [e.index for e in edges]

        targeted_nodes = [e.target for e in edges]
        achieved_tvs = [_graph.es[e]["mono_weight"] + sum_weight for e in eids]
        val_fs = [
            _func_dist(_n_pdb[out_edge.target, :], tv_interval - g_pre)
            for out_edge, g_pre in zip(edges, achieved_tvs)
        ]

        pruned_targeted_nodes = [tn for tn, val_f in zip(targeted_nodes, val_fs) if val_f]
        pruned_achieved_tvs = [tv for tv, val_f in zip(achieved_tvs, val_fs) if val_f]

        for x, y in zip(pruned_targeted_nodes, pruned_achieved_tvs):
            queue.append([x, how_to_get + [cur_node], y])

    return paths


def top_sort_query(start, stop, tv_interval, _graph, _n_pdb):
    """ Retrieve paths using the top. sorted nodes """
    # retrieve top sort first
    # TODO we need to load this from file? Or is it quick enough?
    _top_sort = _graph.topological_sorting()

    dd = defaultdict(lambda: [[], []])

    dd[_top_sort[0]][0] = [0]
    dd[_top_sort[0]][1] = [[_top_sort[0]]]
    for n in _top_sort[0:-1]:

        edges = _graph.es.select(_source=n)
        eids = [e.index for e in edges]
        expand_tvs = [_graph.es[e]["mono_weight"] for e in eids]
        achieved_tvs = [[p_tv + e_tv for e_tv in expand_tvs] for p_tv in dd[n][0]]

        # check if we expand
        val_fs = [
            [_func_dist(_n_pdb[out_edge.target, :], tv_interval - e_tv) for out_edge, e_tv in zip(edges, a_tvs)]
            for a_tvs in achieved_tvs
        ]

        # get all with val_f == true and compact them in our queue
        for p, v_fs, a_tvs in zip(dd[n][1], val_fs, achieved_tvs):
            for e, fs, cur_tv in zip(edges, v_fs, a_tvs):
                if fs:
                    dd[e.target][0].append(cur_tv)
                    dd[e.target][1].append([*p, e.target])

        # save memory
        del dd[n]

    return dd[_top_sort[-1]][1]


def top_sort_attrs_query(start, stop, tv_interval, _graph, _n_pdb):
    """ Retrieve paths using the top. sorted nodes """
    # retrieve top sort first
    # TODO we need to load this from file? Or is it quick enough?
    # This is the fastest imple right now!

    sorted_by_position_attr = []
    s = set([x for x, y in zip(range(_graph.vcount()), _graph.vs.indegree()) if y == 0])
    marked_edges = [False]*_graph.ecount()

    while len(s) != 0:
        t = []  # (isoform_name, iso_pos, pos, n)
        for x in s:
            node_attrs = _graph.vs[x].attributes()

            if "isoform_accession" in node_attrs:
                t1 = node_attrs["isoform_accession"] if node_attrs["isoform_accession"] else node_attrs["accession"]
            else:
                t1 = node_attrs["accession"]

            if "isoform_position" in node_attrs:
                t2 = node_attrs["isoform_position"] if node_attrs["isoform_position"] else float("-inf")
            else:
                t2 = float("-inf")

            t3 = node_attrs["position"] if node_attrs["position"] else float("-inf")
            t.append((t1, t2, t3, x))
        # sorted up down down
        res = sorted(t, key=lambda x: (-len(x[0]), [-ord(c) for c in x[0]], x[1], x[2]))
        n = res[0][3]
        s.remove(n)

        sorted_by_position_attr.append(n)
        for e_out in _graph.vs[n].out_edges():
            marked_edges[e_out.index] = True
            in_out_edges = [x.index for x in _graph.vs[e_out.target].in_edges()]
            if all([marked_edges[x] for x in in_out_edges]):
                s.add(e_out.target)

    _top_sort = sorted_by_position_attr

    dd = defaultdict(lambda: [[], []])

    dd[_top_sort[0]][0] = [0]
    dd[_top_sort[0]][1] = [[_top_sort[0]]]
    for n in _top_sort[0:-1]:

        edges = _graph.es.select(_source=n)
        eids = [e.index for e in edges]
        expand_tvs = [_graph.es[e]["mono_weight"] for e in eids]
        achieved_tvs = [[p_tv + e_tv for e_tv in expand_tvs] for p_tv in dd[n][0]]

        # check if we expand
        val_fs = [
            [_func_dist(_n_pdb[out_edge.target, :], tv_interval - e_tv) for out_edge, e_tv in zip(edges, a_tvs)]
            for a_tvs in achieved_tvs
        ]

        # get all with val_f == true and compact them in our queue
        for p, v_fs, a_tvs in zip(dd[n][1], val_fs, achieved_tvs):
            for e, fs, cur_tv in zip(edges, v_fs, a_tvs):
                if fs:
                    dd[e.target][0].append(cur_tv)
                    dd[e.target][1].append([*p, e.target])

        # save memory
        del dd[n]

    return dd[_top_sort[-1]][1]



def top_sort_attrs_query_limit_variants(start, stop, tv_interval, _graph, _n_pdb, _limit_variants=3):
    """ Retrieve paths using the top. sorted nodes """
    # retrieve top sort first
    # TODO we need to load this from file? Or is it quick enough?
    # This implementation is slower but "SHOULD" reduce the searcch space of p53 considerable!
    # However, the search do take roughly the same amount of time.... (# TODO ...)

    sorted_by_position_attr = []
    s = set([x for x, y in zip(range(_graph.vcount()), _graph.vs.indegree()) if y == 0])
    marked_edges = [False]*_graph.ecount()

    while len(s) != 0:
        t = []  # (isoform_name, iso_pos, pos, n)
        for x in s:
            node_attrs = _graph.vs[x].attributes()

            if "isoform_accession" in node_attrs:
                t1 = node_attrs["isoform_accession"] if node_attrs["isoform_accession"] else node_attrs["accession"]
            else:
                t1 = node_attrs["accession"]

            if "isoform_position" in node_attrs:
                t2 = node_attrs["isoform_position"] if node_attrs["isoform_position"] else float("-inf")
            else:
                t2 = float("-inf")

            t3 = node_attrs["position"] if node_attrs["position"] else float("-inf")
            t.append((t1, t2, t3, x))
        # sorted up down down
        res = sorted(t, key=lambda x: (-len(x[0]), [-ord(c) for c in x[0]], x[1], x[2]))
        n = res[0][3]
        s.remove(n)

        sorted_by_position_attr.append(n)
        for e_out in _graph.vs[n].out_edges():
            marked_edges[e_out.index] = True
            in_out_edges = [x.index for x in _graph.vs[e_out.target].in_edges()]
            if all([marked_edges[x] for x in in_out_edges]):
                s.add(e_out.target)

    _top_sort = sorted_by_position_attr

    dd = defaultdict(lambda: [[], [], []])
    # param dd[key][2] yields number of variants

    dd[_top_sort[0]][0] = [0]
    dd[_top_sort[0]][1] = [[_top_sort[0]]]
    dd[_top_sort[0]][2] = [0]  # Number of variants is 0 at beginning
    for n in _top_sort[0:-1]:

        edges = _graph.es.select(_source=n)
        eids = [e.index for e in edges]
        e_qualifiers = [_resolve_or(x, "VARIANT", min) for x in _graph.es[eids]["qualifiers"]]
        expand_tvs = [_graph.es[e]["mono_weight"] for e in eids]
        achieved_tvs = [[p_tv + e_tv for e_tv in expand_tvs] for p_tv in dd[n][0]]
        acieved_e_qualiefiers = [[var_count + e_q for e_q in e_qualifiers] for var_count in dd[n][2]]

        # check if we expand
        val_fs = [
            [
                # Check for distance
                _func_dist(_n_pdb[out_edge.target, :], tv_interval - e_tv) and
                # Check for variants
                var_count <= _limit_variants
                for out_edge, e_tv, var_count in zip(edges, a_tvs, a_e_vars)
            ]
            for a_tvs, a_e_vars in zip(achieved_tvs, acieved_e_qualiefiers)
        ]

        # get all with val_f == true and compact them in our queue
        for p, v_fs, a_tvs, a_e_vars in zip(dd[n][1], val_fs, achieved_tvs, acieved_e_qualiefiers):
            for e, fs, cur_tv, var_count in zip(edges, v_fs, a_tvs, a_e_vars):
                if fs:
                    dd[e.target][0].append(cur_tv)
                    dd[e.target][1].append([*p, e.target])
                    dd[e.target][2].append(var_count)
                    if var_count == 3:
                        pass

                    if var_count > 3:
                        pass

        # save memory
        del dd[n]

    return dd[_top_sort[-1]][1]

def _resolve_or(fts, feature_type, or_count):
    if fts is None:
        return 0
    count = 0
    for ft in fts:
        if isinstance(ft, Or):
            t = [_resolve_or(or_ft, feature_type, or_count) for or_ft in ft]
            count += or_count(t)
        else:
            if ft.type == feature_type:
                count += 1

    return count

def _func_dist(pdb, s_interval):
    """ function to decide wheather an interval is overlapping or not in pdb. """
    lower_index = np.searchsorted(pdb[:,1], s_interval[0])
    upper_index = np.searchsorted(pdb[:,0], s_interval[1])

    if upper_index > lower_index  or (len(pdb) != lower_index and lower_index == upper_index and pdb[:,0][lower_index] == s_interval[1]):
        # Also check for edge case on the right side of queried s_interval
        return True
    return False