"""Generate Random String"""
from random import choices


def random_string(length):
    """
    The random_string function returns a random string of length n.

    :param length: Determine the length of the random string
    :return: A string of random characters
    :doc-author: Trelent
    """
    abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    return "".join(choices(abc, k=length))
