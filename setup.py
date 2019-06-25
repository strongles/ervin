import subprocess
import setuptools

git_version_binary = subprocess.run(["git", "describe"], capture_output=True).stdout
version = git_version_binary.decode("utf-8").strip()

with open("README.md") as readme:
    long_desc = readme.read()

with open("requirements.txt") as requires:
    require_list = [line.strip() for line in requires.readlines()]

print(require_list)

setuptools.setup(
    name="ervin",
    version=version,
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
    download_url=f"https://github.com/strongles/ervin/archive/{version}.tar.gz",
    packages=setuptools.find_packages(),
    install_requires=require_list,
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
