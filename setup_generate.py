#!/usr/bin/env python3

import argparse
import subprocess


NEWLINE = '\n'
BIG_INDENT = ' ' * 12
SMALL_INDENT = ' ' * 4
SEPARATOR = f",{NEWLINE}{BIG_INDENT}"
SETUP_CONTENT = """import setuptools

with open("README.md") as readme:
    long_desc = str(readme.read())

setuptools.setup(
    name="ervin",
    version="{version}",
    description="ERVin is a collection of tools developed to assist "
                "in discovering ERV sequences within genomic data",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    author="Rob Williams",
    author_email="robccwilliams@gmail.com",
    entry_points={{
        "console_scripts": ["ervin = ervin.ervin:run"]
    }},
    licence="GPL-3.0-or-later",
    url="https://github.com/strongles/ervin",
    download_url="https://github.com/strongles/ervin/archive/{version}.tar.gz",
    packages=setuptools.find_packages(),
    install_requires={require_list},
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
"""


def get_current_git_tag():
    git_version_binary = subprocess.run(["git", "describe", "--abbrev=0"], capture_output=True).stdout
    return git_version_binary.decode("utf-8").strip()


def get_requirements():
    with open("requirements.txt") as requires:
        require_list = [f'"{line.strip()}"' for line in requires.readlines()]
    return f"[{NEWLINE}{BIG_INDENT}{SEPARATOR.join(require_list)}{NEWLINE}{SMALL_INDENT}]"


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--minor-version", action="store_true", dest="minor_version")
    parser.add_argument("--mid-version", action="store_true", dest="mid_version")
    parser.add_argument("--major-version", action="store_true", dest="major_version")
    return parser.parse_args()


def build_version_number(arguments):
    current_git_tag = get_current_git_tag()
    tag_components = [int(component) for component in current_git_tag.split(".")]
    if arguments.major_version:
        tag_components[0] += 1
    if arguments.mid_version:
        tag_components[1] += 1
    if arguments.minor_version:
        tag_components[2] += 1
    return ".".join([str(component) for component in tag_components])


if __name__ == "__main__":
    args = get_args()
    new_version_number = build_version_number(args)
    requirements = get_requirements()
    with open("setup.py", 'w') as setup_file:
        setup_file.write(SETUP_CONTENT.format(version=new_version_number, require_list=requirements))
