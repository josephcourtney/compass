#!/usr/bin/env python
from distutils.core import setup, Extension

setup(
    name = "compass",
    version = "0.1.1",
    description = 'COMPASS structure determination package',
    author = 'Rienstra Group, UIUC',
    packages = ['compass'],
    py_modules = [
        'compass.ros', 
        'compass.mod'
    ],
    ext_modules = [
        Extension(
            'compass.util',
            sources = ['compass/util.c']
        )
    ]
)
