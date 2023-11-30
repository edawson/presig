import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
        name="presig",
        version="0.1.0",
        author="Eric T. Dawson",
        author_email="eric@betulalabs.com",
        description="presig: a package for preparing MAF files for mutational signature extraction",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/edawson/presig",
        packages=setuptools.find_packages(exclude=[]),
        install_requires=["pyfaidx>=0.5.5.2"],
        entry_points = {
            'console_scripts': [
                'spectre = presig.presig:main',
                ],
            },
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            ],
        python_requires='>=3.6',
        )
