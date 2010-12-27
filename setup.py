from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(name='rnaseq',
      version='0.1',
      description='RNASeq analysis tool',
      author='Fred Ross',
      author_email='madhadron@gmail.com',
      cmdclass={'build_ext': build_ext},
      package=['rnaseq'],
      ext_modules=[Extension("rnaseq.model",["rnaseq/model.pyx"])],
      scripts=['bin/samfiles_to_sqlite.py', 'bin/find_subproblems.py', 
               'bin/inference.py', 'bin/analyze_rnaseq.py'],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics']
      )
