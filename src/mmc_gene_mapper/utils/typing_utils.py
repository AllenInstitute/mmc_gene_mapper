"""
Define function to "simulate" strong typing
"""


def check_many_types(name_arr, arg_arr, type_arr):
    """
    Check that a list of arguments are of the expected typ.
    Raise a ValueError if any are not.

    Parameters
    ----------
    name_arr:
        list of strings. The names of the arguments
        (for use in error message)
    arg_arr:
        list of the arguments being typed
    type_arr:
        list of the types the args in arg_arr
        are expected to be
    """
    msg = ""
    for name, arg, expected_type in zip(name_arr, arg_arr, type_arr):
        msg += check_type(
            arg_name=name,
            arg=arg,
            expected_type=expected_type
        )
    if len(msg) > 0:
        raise ValueError(msg)


def check_type(arg_name, arg, expected_type):
    """
    Parameters
    ----------
    arg_name:
        str; the name of the parameter being typed
    arg:
        the value of the paramter being typed
    expected_type:
        the type arg is expected to be

    Returns
    -------
    A str; the error message (if any) indicating why the typing
    is wrong
    """
    msg = ""
    if not isinstance(arg, expected_type):
        msg = (
            f"{arg_name} must be of type {expected_type}; "
            f"you gave {arg} of type {type(arg)}\n"
        )
    return msg
