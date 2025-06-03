import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale_nanopore", # Replace with your own username
    version="1.0.1",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data - Nanopore application",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/apscale_nanopore/",
    packages=setuptools.find_packages(),
    license='MIT',
    install_requires=[
                    'Bio >= 1.7.1',
                    'biopython >= 1.85',
                    'joblib >= 1.4.2',
                    'ete3 >= 3.1.3',
                    'numpy >= 1.26.4',
                    'pandas >= 2.2.3',
                    'pyarrow >= 19.0.0',
                    'xmltodict >= 0.14.2',
                    'plotly >= 6.1.1',
                    'tqdm >= 4.67.1',
                ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    entry_points={
        "console_scripts": [
            "apscale_nanopore = apscale_nanopore.__main__:main",
        ]
    },
)

