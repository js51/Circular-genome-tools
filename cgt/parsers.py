def newick_to_json(
    newick_string,
    namestring="id",
    lengthstring="branch_length",
    childrenstring="children",
    generate_names=False,
):
    """Parses a newick string (example "((A,B),C)") into JSON format accepted by NetworkX

    Arguments:
            newick_string: The newick string to convert to JSON. Names, branch lengths and named
                    internal nodes accepted. If no names are given, an ID is generated for the name if generate_names=True.
            namestring: The label for node name to use in the JSON (default "id").
            lengthstring: The label for branch length to use in the JSON (default "branch_length").
            childrenstring: The label for children arrays to use in the JSON (default "children").
            generate_names: Whether or not to generate names if none are given
    """
    newick_string = newick_string.strip(";")
    node = dict()
    info = newick_string[newick_string.rfind(")") + 1:]
    info = info.split(":")
    name = info[0]
    node[namestring] = name
    length = ""
    if len(info) > 1:
        length = info[1]
        node[lengthstring] = float(length)
    children_string = newick_string[0 : newick_string.rfind(")") + 1]
    if children_string != "":
        if children_string[-1] == ")":
            children_string = children_string[1:-1]
        children_string = __split_into_children(children_string)
        child_nodes_JSON = []
        generated_name_parts = []
        for child_string in children_string:
            child_JSON = newick_to_json(
                child_string, namestring, lengthstring, childrenstring, generate_names
            )
            if (not name) and generate_names:
                generated_name_parts.append(child_JSON[namestring])
            child_nodes_JSON.append(child_JSON)
        if generate_names:
            node[namestring] = "|".join(generated_name_parts)
        node[childrenstring] = child_nodes_JSON
    return node


def json_to_newick(json_dict):
    """Converts a JSON tree to a newick string

    Arguments:
            json_dict: The JSON tree to convert to newick.
    """
    return __json_to_newick(json_dict) + ";"


def __json_to_newick(
    json_dict, namestring="id", lengthstring="branch_length", childrenstring="children"
):
    if not isinstance(json_dict, list):
        string = ""
        if childrenstring in json_dict:
            string += (
                "("
                + __json_to_newick(json_dict[childrenstring])
                + ")"
                + json_dict[namestring]
            )
        else:
            string += json_dict[namestring]
        if lengthstring in json_dict:
            string += ":" + str(json_dict[lengthstring])
        return string
    else:
        string = ""
        for dict in json_dict:
            string += __json_to_newick(dict) + ","
        return string[0:-1]


def __split_into_children(children_string):
    """Helper function for splitting a newick string into siblings"""
    string = children_string
    opened = 0
    closed = 0
    pieces = []
    last_stop = 0
    for i, s in enumerate(string):
        if s == "(":
            opened += 1
        if s == ")":
            closed += 1
        if s == "," and opened == closed:
            pieces.append(string[last_stop:i])
            last_stop = i + 1
    pieces.append(string[last_stop:])
    return pieces
