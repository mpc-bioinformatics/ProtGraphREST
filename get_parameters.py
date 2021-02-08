from sanic import exceptions

def get_path(request):

    paths = []

    if "path" in request.args:
        path = request.args["path"][0]
        # Parse the integer list
        try:
            path_ints = list(map(int, path.split(",")))
            paths.append(path_ints)
        except Exception:
            raise exceptions.ServerError("'path' can only consist of ',' and '[0-9]'. Actual: '{}'".format(path), status_code=400)


    if request.method == "POST":
        if request.json is not None and "path" in request.json:
            path = request.json["path"]
            if not all(isinstance(x, int) for x in path):
                raise exceptions.ServerError("'path' can only be a list of integer values! Actual: '{}'".format(path), status_code=400)
            paths.append(path)


    if len(paths) == 0:
        pass  # Return error no input given!







    request.args



    print("yeah")


