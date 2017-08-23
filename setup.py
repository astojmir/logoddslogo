#!/usr/bin/env python

import sys

from distutils.core import setup
from distutils.core import Extension
from distutils.command.build import build
from distutils.command.install_data import install_data

# Supress warning that distutils generates for the install_requires option
import warnings
warnings.simplefilter('ignore', UserWarning, lineno=236)

# check dependencies
if not hasattr(sys, 'version_info') or sys.version_info < (2,7,0,'final'):
    raise SystemExit("Dependency error: LogOddsLogo requires Python 2.7 or later.")


from logoddslogolib import __version__

EXTENSIONS = [
    Extension('logoddslogolib._logodds',
              ['logodds/info1.c',
               'logodds/bildn.c',
               'logodds/bildp.c',
               'logodds/logodds.c',
               'logodds/weighcounts.c',
               'logodds/schneider.c',
               ],
              libraries=['m'],
              ),
    ]

def main() :
    long_description = open("README.txt").read()


    setup(
        name             =  "logoddslogo",
        version          =  __version__,
        description      = "LogOddsLogo : Log Odds Sequence Logos",
        long_description  = long_description,
        maintainer       = "Yi-Kuo Yu",
        maintainer_email = "yyu@ncbi.nlm.nih.gov",
        url              = "http://www.ncbi.nlm.nih.gov/CBBresearch/Yu/logoddslogo/",

        download_url     = 'ftp://ftp.ncbi.nlm.nih.gov/pub/qmbp/logoddslogo/logoddslogo-%s.tar.gz' % __version__ ,
        classifiers      = [
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'License :: Public Domain',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Topic :: Software Development :: Libraries',
            'Topic :: Software Development :: Libraries :: Python Modules',
            ],

        scripts         = ['logoddslogo', 'transformseq'],
        # ext_modules=EXTENSIONS,
        packages  = ['corebio',
                     'corebio.db',
                     'corebio.secstruc',
                     'corebio.seq_io',
                     'corebio.seq_io._nexus',
                     'corebio.ssearch_io',
                     'corebio.utils',
                     'weblogolib',
                     'logoddslogolib',
                     ],
        package_data={'weblogolib': ['htdocs/*.*',
                                     'htdocs/img/*.*',
                                     'template.eps'],
                      'logoddslogolib': ['htdocs/*.*',
                                         'htdocs/examples/*.*',
                                         'htdocs/P/*.*',
                                         'htdocs/PF/*.*',
                                         'htdocs/qmbppics/*.*',
                                         'htdocs/apidocs/*.*',
                                         'template.eps'],
                      'corebio': ['data/*.*']
                      },
        requires=['numpy'],
    )


if __name__ == '__main__' :
    main()
