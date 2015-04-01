#!/usr/bin/env python

'''
setup.py file for groupfilter
'''

from distutils.core import setup

version = open('VERSION').readline().strip()

setup(
    name = 'groupfilter',
    version = version,
    description      = '''A software tool for efficient filtering of Morpheus search engine results.''',
    long_description = (''.join(open('README').readlines())),
    author           = 'Mark Ivanov & Lev Levitsky',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'http://hg.theorchromo.ru/groupfilter',
    requires         = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Topic :: Education',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Chemistry',
                        'Topic :: Scientific/Engineering :: Physics',
                        'Topic :: Software Development :: Libraries'],
    license          = 'License :: OSI Approved :: Apache Software License',
    packages = ['groupfilter'],
    package_dir = {'groupfilter': 'groupfilter'},
    scripts=['groupfilter/groupfilter.py']
    )