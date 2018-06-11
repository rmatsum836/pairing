import os
from pkg_resources import resource_filename


def get_fn(name):
    """Get the full path to one of the reference files shipped for utils.
    In the source distribution, these files are in ``pairing/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the reference/ folder).

    Returns
    -------
    fn : str
        Full path to file to load, including name
    """

    fn = resource_filename('pairing', os.path.join('data', name))
    if not os.path.exists(fn):
        raise IOError('Sorry! {} does not exists.'.format(fn))
    return fn
