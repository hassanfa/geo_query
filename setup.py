from setuptools import setup, find_packages

setup(
    name='geofetch',
    version='1.0.0',  # Update this as needed
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

