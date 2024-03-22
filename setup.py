#!/usr/bin/env python

from setuptools import setup, find_packages
import os

# https://github.com/ninjaaron/fast-entry_points
import fastentrypoints

with open(
    os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.md'),
    encoding='utf-8',
) as readme_md:
    long_description = readme_md.read()

extras_require = {
    'develop': ['bumpversion', 'black', 'pyflakes'],
    'plotting': ['matplotlib'],
    'test': ['pytest~=5.0', 'pytest-cov>=2.5.1', 'coverage>=4.0', 'pytest-mock'],
    'docs': ['sphinx>=4.0.0', 'sphinx_rtd_theme'],
}
extras_require['complete'] = sorted(set(sum(extras_require.values(), [])))


setup(
    name='pyssrl',
    version='0.0.1',
    package_dir={"": "src"},
    packages=find_packages(where="src", exclude=["tests"]),
    include_package_data=True,
    description='Python analysis for SSRL test beam',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gitlab.cern.ch/scipp/collinearw',
    author='SCIPP',
    author_email='scipp@cern.ch',
    license='',
    keywords='SSRL LGAD analysis',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    python_requires=">=3.8",
    install_requires=[
        'click',
        'formulate',
        'tabulate',
        'numexpr',
        'numpy',
        'uproot',
        'awkward',
        'tqdm',
        'pandas',
        'rich',
        'numba',
        'lazy_loader',
    ],
    extras_require=extras_require,
    dependency_links=[],
    entry_points={'console_scripts': ['pyssrl=pyssrl.cli:pyssrl']},
)

