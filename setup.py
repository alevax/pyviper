import setuptools
  
with open("README.md", "r") as fh:
    description = fh.read()
  
setuptools.setup(
    name="Pyther",
    version="1.0.0",
    author="Alexander LE Wang & Zizhao Lin",
    author_email="aw3436@cumc.columbia.edu",
    packages=["Pyther"],
    description="A package to load ARACNe networks and to run VIPER and NaRnEA in Python",
    long_description="This Scanpy-compatible package provides the Interactome class to easily allow the user to both load and interact with ARACNe networks. Once these networks have been loaded, the user has the option to filter regulators, filter targets, and prune the network. The user can then give their Interactome (or multiple Interactomes) and a gene expression signature to Pyther function to either run VIPER or NaRnEA to infer protein activity. This package also provides the option to compute pathway enrichment with Interactomes built from sets of pathways from MSigDB.",
    long_description_content_type="text/markdown",
    url="https://github.com/alevax/pyther",
    license='MIT',
    python_requires='>=3.8',
    install_requires=[scanpy, scipy, os, shutil, joblib, pathlib]
)
