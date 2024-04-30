import setuptools

with open("README.md", "r", encoding = "utf-8") as fh:
    description = fh.read()

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setuptools.setup(
    name="viper-in-python",
    version="1.0.8",
    author="Alexander L.E. Wang & Zizhao Lin & Luca Zanella",
    author_email="aw3436@cumc.columbia.edu",
    packages=setuptools.find_packages(),
    package_data={'pyviper': ['data/*']},  # Include data directory
    #package_data={'': ['data/*']}, # Include data directory
    include_package_data=True,
    description="A package for VIPER-based Protein Activity analysis of transcriptomic data in Python",
    long_description=long_description,
    #long_description="This Scanpy-compatible package provides the Interactome class to easily allow the user to both load and interact with gene regulatory networks. Once these networks have been loaded, the user has the option to filter regulators, filter targets, and prune the network. The user can then give their Interactome (or multiple Interactomes) and a gene expression signature to the viper function to run VIPER (aREA or NaRnEA enrichment methods) to infer protein activity. This package also provides the option to compute pathway enrichment with Interactomes built from sets of pathways from MSigDB and several modules for data processing.",
    long_description_content_type="text/markdown",
    url="https://github.com/alevax/pyviper",
    project_urls = {
    	'Documentation': 'https://alevax.github.io/pyviper/',
    	'Source': 'https://github.com/alevax/pyviper',
    },
    license='MIT',
    python_requires='>=3.8',
    install_requires=[
        "scipy",
        "tqdm",
        "scanpy",
        "anndata",
        "pandas>=1.3.0",
        "pandas<2.0",
        "numpy",
        "joblib",
        "statsmodels"
    ]
)
