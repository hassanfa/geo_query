from setuptools import setup, find_packages

# Read the version from version.txt
with open('version.txt') as version_file:
    version = version_file.read().strip()

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

