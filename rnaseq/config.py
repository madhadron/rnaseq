import os
from ConfigParser import SafeConfigParser, NoOptionError
from string import lstrip, rstrip

class ConfigurationError(Exception):
    pass

def to_absolute(filename):
    """If *filename* refers to a file that exists, return its absolute path.

    If *filename* does not exist, raises a ConfigurationError.
    """
    abs_filename = os.path.abspath(filename)
    if os.path.exists(abs_filename):
        return abs_filename
    else:
        raise ConfigurationError("File %s in configuration doesn't exist." % filename)

def load_configuration(config_file):
    cp = SafeConfigParser({})
    cp.read(config_file)
    r = {}
    for s in cp.sections():
        try:
            files = [lstrip(rstrip(f)) for f in 
                     cp.get(s, 'fastqfiles').split(',')]
            abs_files = [to_absolute(f) for f in files]
            r[s] = {'control': cp.getboolean(s, 'control'),
                    'fastqfiles': abs_files}
        except NoOptionError, n:
            raise ConfigurationError("Group '%s' in %s is missing '%s' option" % (n.section,
                                                                                  config_file,
                                                                                  n.option))
    return r
