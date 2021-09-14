#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="conflict",
    version="0.5",
    description="stochastic, gridded conflict model",
    author="Scott D. Peckham",
    author_email="speckham@air.org",
    license="MIT",
    url="https://github.com/peckhams/stochastic_conflict_model",
    include_package_data=True,
    packages=find_packages("."),
    install_requires=["numpy", "gdal", "pandas", "matplotlib", "pip"],
)
