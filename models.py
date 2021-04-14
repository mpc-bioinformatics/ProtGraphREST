from typing import List, Literal, Optional

from pydantic import BaseModel, validator


class Path(BaseModel):
    """ BaseClass, which is used to parse one or multiple paths into list of ints """

    # Parse it first to List of ints, then strings (comma or -> seperated)
    path: Optional[List[int]] = []
    # Parse a list of strings (comma or -> seperated), then strings (; sperated (comma or -> seperated))
    paths: Optional[List[List[int]]] = []
    # Set the return type, we always can return json or text here (text -> DEFAULT)
    returns: Optional[Literal["text", "json"]] = "text"

    @validator("path", pre=True)
    def convert_str_to_list(cls, v):  # pylint: disable=E0213
        """ Allow string paths seperated by '->' or ',' """
        if type(v) == str:
            return v.replace("->", ",").split(",")
        return v

    @validator("paths", pre=True)
    def convert_str_to_list_of_lists(cls, v):  # pylint: disable=E0213
        """ Allow string paths seperated as above and by ';' """
        if type(v) == str:
            return [x.replace("->", ",").split(",") for x in v.split(";")]
        elif type(v) == list:
            for idx, x in enumerate(v):
                if type(x) == str:
                    v[idx] = x.replace("->", ",").split(",")
        return v


class MonoWeigthQuery(BaseModel):
    """ BaseClass, which is used to parse one or multiple parameters for weight queries """

    # Parse the unit, either ppm or Da, which is used to calculate the interval
    unit: Literal["ppm", "Da"]

    # Parse the mass tolerance around the mono_weight
    mass_tolerance: int

    # The weight, which we want to query
    mono_weight: float

    # The time we want to retrieve the result (in seconds). If it takes longer than timeout, we abort the search
    timeout: Optional[float] = 10000

    # Set the number of intervals to be used
    k: Optional[int] = 10
