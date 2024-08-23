from setuptools import setup, find_packages

setup(
    name='cdm_utils',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # List your package dependencies here, e.g.,
        # 'numpy>=1.18.0',
        # 'pandas>=1.0.0',
    ],
    entry_points={
        'console_scripts': [
            # Define command-line scripts here if any
        ],
    },
    test_suite='tests',
)

