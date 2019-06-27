import setuptools

with open("README.md") as readme:
    long_desc = str(readme.read())

setuptools.setup(
    name="ervin",
    version="0.0.5",
    description="ERVin is a collection of tools developed to assist "
                "in discovering ERV sequences within genomic data",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    author="Rob Williams",
    author_email="robccwilliams@gmail.com",
    entry_points={
        "console_scripts": ["ervin = ervin.ervin:run"]
    },
    licence="GPL-3.0-or-later",
    url="https://github.com/strongles/ervin",
    download_url="https://github.com/strongles/ervin/archive/0.0.5.tar.gz",
    packages=setuptools.find_packages(),
    install_requires=[
            "biopython==1.73",
            "numpy==1.16.1",
            "progressbar2==3.39.3",
            "nose==1.3.7",
            "mock==3.0.5",
            "flake8==3.7.7",
            "ftputil==3.4"
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Environment :: MacOS X",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
