"""A setuptools based setup module.
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='karenina',
    version='0.0.1',
    description="This script simulates microbiome " +
    "change over time using Ornstein-Uhlenbeck (OU) models.  These are " +
    "similar to Brownian motion models, with the exception that they " +
    "include reversion to a mean. Output is a tab-delimited data table " +
    "and figures.",
    url='https://github.com/zaneveld/karenina',
    author='Jesse Zaneveld',
    author_email='zaneveld@gmail.com',
    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Qiime2 Users',
        'Topic :: Microbiology :: Visualization and Modeling Tools',
        'License :: GPL',
        #'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.7',
        #'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        #'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='karenina ornstein uhlenbeck brownian motion',
    packages=find_packages(exclude=['docs', 'tests', 'q2-karenina']),
    install_requires=['scipy'],
    package_data={
        "karenina.data":['add_x_mu_high.tsv', 'all_perturbations.tsv',
                  'double_xyz_delta.tsv', 'double_z_delta.tsv',
                  'README.txt', 'set_x_lambda_high.tsv',
                  'set_x_lambda_medium.tsv', 'set_x_lambda_small.tsv',
                  'set_x_mu_high.tsv', 'set_xyz_mu_high.tsv',
                  'set_xyz_mu_low.tsv', 'set_y_lambda_high.tsv',
                  'set_yz_lambda_high.tsv', 'set_z_lambda_zero.tsv',
                  'set_x_mu_low.tsv', 'set_xyz_lambda_low.tsv',
                  'set_xyz_lambda_zero.tsv', 'set_xyz_lambda_zero.tsv',
                  'set_y_lambda_medium.tsv', 'set_yz_lambda_medium.tsv']},

    entry_points={
        'console_scripts': [
            'spatial_ornstein_uhlenbeck=karenina.spatial_ornstein_uhlenbeck:main',
            'fit_timeseries=karenina.fit_timeseries:main',
            'benchmark = karenina.benchmark:main',
            'karenina_visualization = karenina.visualization:main'
        ],
    },
)
