groupfilter - a software tool for efficient filtering of Morpheus search engine results
---------------------------------------------------------------------------------------

**groupfilter** is a module based on the
**pyteomics** (https://pyteomics.readthedocs.io/) package that
increases the number of identifications above FDR threshold for Morpheus search engine results.

Dependencies
------------

- pyteomics (>= 3.0)
- numpy

How to install
--------------

 - best way to install **groupfilter** is with ``pip``: ``pip install numpy
   pyteomics groupfilter``.
 - alternatively, you can download source code for each library and run
   ``python setup.py install`` for each of them.
 - There may be other ways to install ``numpy``, depending on your platform.

How to use
--------------

``groupfilter.py file1.tsv [file2.tsv ...] db.fasta``

positional arguments:
  files       list of Morpheus PSM files

  db.fasta       path to FASTA file

Links
-----

- PyPI: https://pypi.python.org/pypi/groupfilter
- BitBucket repo & issue tracker: https://github.com/levitsky/groupfilter
- Mailing list: pyteomics@googlegroups.com
