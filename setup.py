#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='seisvo',
    version='1.0.0',
    description='A python tool for seismovolcanic processes',
    author='Ivan Melchor',
    author_email='ifmelchor@unrn.edu.ar',
    url='https://github.com/ifmelchor/seisvo',
    packages=find_packages(),
    install_requires=[
        'obspy>=1.3.0',
        'geopy>=2.3',
        'utm>=0.7.0',
        'h5py>=3.7',
        'juliacall>=0.9.11',
        'py-notifier>=0.5.0',
        'pyqtgraph>=0.12.4',
        'simplekml>=1.3.6',
        'geopy>=0.5.0',
        'scikit-learn>=1.1.2',
        'pandas>=1.4.3',
        'matplotlib>=3.8',
    ],
    python_requires=">=3.8",
)
