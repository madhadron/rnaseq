from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(name='rnaseq',
      version='0.1',
      description='RNASeq analysis pipeline',
      author='Fred Ross',
      author_email='madhadron@gmail.com',
      cmdclass={'build_ext': build_ext},
      modules=['bag.py'],
      ext_modules=[Extension("rnaseq",["rnaseq.pyx"])],
      scripts=['samfiles_to_sqlite.py', 'find_subproblems.py', 
               'inference.py', 'analyze_rnaseq.py'],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics']
      )
