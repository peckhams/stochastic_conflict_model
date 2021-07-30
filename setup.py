#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="conflict",
    version="0.5",
    description="stochastic, localized conflict model",
    author="Scott D. Peckham",
    author_email="scott.peckham@kimetrica.com",
    license="MIT",
    url="https://github.com/peckhams/conflict",
    include_package_data=True,
    packages=find_packages("."),
    install_requires=["numpy", "gdal", "osr", "pandas"],
)
