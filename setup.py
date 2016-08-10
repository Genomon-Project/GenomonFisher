#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_fisher',
    version='0.2.0',
    description='Python tools to identify somatic mutations.',
    author='Ken-ichi Chiba',
    author_email='kchiba@hgc.jp',
    url='https://github.com/Genomon-Project/GenomonFisher',
    package_dir = {'': 'lib'},
    packages=['fisher','bamfilter'],
    scripts=['fisher','bamfilter'],
    license='GPL-3'
)
