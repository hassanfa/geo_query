from setuptools import setup, find_packages
import os


package_root = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(package_root, "geo_query/version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

setup(
    name='geofetch',
    version=version,
    description='A simplistic command-line tool to fetch GEO data based on user input.',
    author='Hassan Foroughi',
    packages=find_packages(),
    install_requires=open('requirements.txt').read().splitlines(),  # Reads from requirements.txt
    entry_points={
        'console_scripts': [
            'geofetch=geo_query.geo_query:cli',  # Adjusted to match your structure
        ],
    },
    python_requires='>=3.6',
)

