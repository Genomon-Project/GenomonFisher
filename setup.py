#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_fisher',
    version='0.1.0',
    description='Python tools to identify somatic mutations.',
    author='Eigo Shimizu',
    author_email='eigos@hgc.jp',
    url='https://github.com/Genomon-Project/GenomonFisher',
    package_dir = {'': 'lib'},
    packages=['fisher'],
    scripts=['fisher'],
    license='GPL-3'
)
