# setup.py file for the wrpoly package
# Authors: Danielle Alverson, Eric Fonseca, Kausturi Parui, Talianna Ulloa, and Bonnie Stolt
from setuptools import setup

setup(
    name="wrpoly",
    version="0.1.0",
    author="Danielle Alverson, Eric Fonseca, Kausturi Parui, Talianna Ulloa, and Bonnie Stolt",
    author_email="ericfonseca@ufl.edu", 
    description="A package for calculating the polyhedra of waldley-roth complexes",
    license='BSD 2-clause',
    packages=['wrpoly'],
    install_requires=['numpy', 'scipy', 'matplotlib', 'pymatgen', 'pandas'],
    package_data={'wrpoly': ['data/*.csv']},

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
