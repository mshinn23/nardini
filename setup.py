from setuptools import setup


description = ' '.join("""Nardini is a Python package for IDP sequences
that generates similar sequences derived from a statistical analysis of
residue-specific asymmetries across short blob sizes, calculated across
bootstrapped sequences.""".split())


setup(
    name='nardini',
    description=description,
    version='1.1',
    author='Min Kyung Shinn, Kiersten Ruff',
    author_email='mshinn@wustl.edu, kiersten.ruff@wustl.edu',
    maintainer='Jared Lalmansingh',
    maintainer_email='jared.lalmansingh@wustl.edu',
    license='GPL v2',
    url='https://github.com/mshinn23/nardini',
    packages=['nardini'],
    scripts=['bin/nardini'],
    python_requires='>=3.5',
    install_requires=[
        'BioPython>=1.75',
        'matplotlib>=3.3',
        'numpy>=1.19',
        'scipy<1.7',
        'pandas>=1.2',
        'tabulate>=0.8'
    ],
    include_package_data=True,
    zip_safe=True
)
