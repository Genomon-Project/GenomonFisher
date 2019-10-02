#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='genomon_fisher',
    version='0.4.3',
    description='Python tools to identify somatic mutations.',
    author='Ken-ichi Chiba',
    author_email='kchiba@hgc.jp',
    url='https://github.com/Genomon-Project/GenomonFisher',
    # package_dir = {'': 'lib'},
    packages = find_packages(exclude = ['tests']),
    test_suite = 'unit_tests.suite',
    license='GPL-3'
)

