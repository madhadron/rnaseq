"""
Configuration tests.

>>> from rnaseq.config import *

>>> to_absolute("/bin/sh")
'/bin/sh'

>>> to_absolute("test_data")[-9:]
'test_data'

>>> to_absolute("/dev/null/boris")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "rnaseq/config.py", line 18, in to_absolute
    raise ConfigurationError("File %s in configuration doesn't exist." % filename)
ConfigurationError: File /dev/null/boris in configuration doesn't exist.





"""

if __name__ == '__main__':
    import doctest
    doctest.testmod()
