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
        'obspy>=1.3.0'
    ],
    python_requires=">=3.8",
)
