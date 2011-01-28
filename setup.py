from distutils.core import setup
from distutils.extension import Extension

numpy_include_dirs = ['/usr/include/python2.6',
                      '/usr/lib64/python2.6/site-packages/numpy/core/include']

try:
    from Cython.Distutils import build_ext
except:
    print "You don't seem to have Cython installed. Please get a"
    print "copy from www.cython.org and install it"
    sys.exit(1)

setup(name='rnaseq',
      version='0.1',
      description='RNASeq analysis tool',
      author='Fred Ross',
      author_email='madhadron@gmail.com',
      cmdclass={'build_ext': build_ext},
      packages=['rnaseq'],
      ext_modules=[Extension("rnaseq.model",["rnaseq/model.pyx"],
                             include_dirs = numpy_include_dirs + ['.'],
                             extra_compile_args = ['-O3', '-Wall'],
                             extra_link_args = ['-g'])],
      scripts=['bin/samfiles_to_sqlite.py', 'bin/find_subproblems.py', 
               'bin/inference.py', 'bin/analyze_rnaseq.py',
               'bin/simple_inference.py', 'bin/inference_subproblem.py',
               'bin/workflow.py'],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics']
      )
