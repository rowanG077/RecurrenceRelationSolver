#!/usr/bin/env python3
# coding=utf-8

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='RecurrenceRelationSolver',
    version='0.1.0',
    description='Solves recurrence relation into a closed form',
    long_description=readme,
    author='Rowan Goemans',
    author_email='rgoemans@science.ru.nl',
    url='https://github.com/rowanG077/RecurrenceRelationSolver',
    license=license,
    packages=find_packages(exclude=('tests')),
    entry_points={
        'console_scripts': [
        'recurrenceSolver=RecurrenceRelationSolver.RecurrenceRelationSolver'
        ]
    },
)
