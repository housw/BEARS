#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['click>=6.0,<=6.7',
                'click-default-group>=1.2.1',
                'pandas>=0.24.0',
                'scipy>=1.3.0',
                'scikit-learn>=0.21.2',
                'biopython>=1.68',
                'colorama>=0.4.1',
                'umap-learn>=0.3.9',
                'hdbscan>=0.8.22',
                'freqgen>=0.1.0',
                'MulticoreTSNE>=0.1'
                ]

setup_requirements = ['pytest-runner',
                      'cython>=0.29.12',
                      'numpy>=1.16.4'
                      ]

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
    long_description=readme + '\n\n',
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

