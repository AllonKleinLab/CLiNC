import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="clinc", # Replace with your own username
    version="0.1.1",
    author="Caleb Weinreb",
    author_email="calebsw@gmail.com",
    description="Cell Lineage from Normalized Covariance ",
    long_description='Cell-Lineage-from-Normalized-Covariance (CLiNC) is a method to reconstruct developmental hierarchies from clonal barcoding data.  Briefly, the model underlying CLiNC assumes that all barcodes are deposited as a synchronous moment in differentiation and that differentiation events are not directly coupled to cell division (as in asymmetric division)',
    long_description_content_type="text/markdown",
    url="https://github.com/AllonKleinLab/Cell-Lineage-from-Normalized-Covariance",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "ete3",
        "numpy",
        "scipy",
        "matplotlib",
        "ipykernel",
        "jupyter",
        "fastcluster",
        "SetCoverPy"
    ],
    python_requires='>=3.6',
) 
