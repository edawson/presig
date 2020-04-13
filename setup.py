import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="edawson", # Replace with your own username
    version="0.0.1",
    author="Eric T. Dawson",
    author_email="eric@erictdawson.com",
    description="presig: a package for preparing MAF files for signature extraction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/edawson/presig",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
