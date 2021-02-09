import falcon
from pydantic import BaseModel, ValidationError

from prot_graph_exception import ProtGraphException


def load_model(model: BaseModel, *params: dict):
    """
    Loads model information, depending on the added params.
    Returns as many objects, as has been added to this function
    """
    results = []
    for x in params:
        try:
            # Parse only the added parameters
            path_obj = model(**x)
            results.append(path_obj)
        except ValidationError as ve:
            raise ProtGraphException(falcon.HTTP_400, ve.json())

    return results[0] if len(results) == 1 else results
