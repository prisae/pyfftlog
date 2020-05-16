# -*- coding: utf-8 -*-
import os
import re
from setuptools import setup

# Get README and remove badges.
readme = open('README.rst').read()
readme = re.sub('----.*marker', '----', readme, flags=re.DOTALL)

setup(
    name='pyfftlog',
    description='Logarithmic Fast Fourier Transform',
    long_description=readme,
    author='Dieter Werthmüller',
    author_email='dieter@werthmuller.org',
    url='https://github.com/prisae/pyfftlog',
    license='CC0-1.0',
    packages=['pyfftlog', ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    ],
    install_requires=[
        'scipy',
    ],
    use_scm_version={
        'root': '.',
        'relative_to': __file__,
        'write_to': os.path.join('pyfftlog', 'version.py'),
    },
    setup_requires=['setuptools_scm'],
)
