#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['click>=6.0,<=6.7',
                'click-default-group>=1.2.1',
                'numpy>=1.16.4',
                'pandas>=0.24.0',
                'scikit-learn>=0.21.2',
                'scipy>=1.3.0',
                'colorama>=0.4.1',
                'biopython>=1.68',
                'pymer>=0.3.0',
                'hdbscan>=0.8.22',
                'freqgen>=0.1.0',
                'bokeh>=1.2.0',
                'Jinja2>=2.10.1',
                'MarkupSafe>=1.1.1',
                'packaging>=19.0',
                'pyparsing>=2.4.0',
                'Pillow>=6.1.0',
                'python-dateutil>=2.8.0',
                'PyYAML>=5.1.1',
                'tornado>=6.0.3',
                'dit>=1.2.3',
                'boltons>=19.1.0',
                'bz2file>=0.98',
                'cffi>=1.12.3',
                'contextlib2>=0.5.5',
                'Cython>=0.29.12',
                'debtcollector>=1.21.0',
                'decorator>=4.4.0',
                'h5py>=2.9.0',
                'joblib>=0.13.2',
                'llvmlite>=0.29.0',
                'MulticoreTSNE>=0.1',
                'networkx>=2.3',
                'numba>=0.44.1',
                'pbr>=5.4.0',
                'prettytable>=0.7.2',
                'pycparser>=2.19',
                'pymer>=0.3.3',
                'pytz>=2019.1',
                'screed>=1.0',
                'six>=1.12.0',
                'umap-learn>=0.3.9',
                'wrapt>=1.11.2'
                ]


setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Shengwei Hou",
    author_email='housw2010@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Binning metagenomic contigs via length-dependent iterative clustering and integration",
    entry_points={
        'console_scripts': [
            'blendit=blendit.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='blendit',
    name='blendit',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/housw/blendit',
    version='0.1.0',
    zip_safe=False,
)
