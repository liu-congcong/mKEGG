import os
from setuptools import setup, find_packages


def main():
    setup(
        name = 'mkegg',
        version = '1.0.0',
        url = 'https://github.com/liu-congcong/mkegg/',
        author = 'Liucongcong',
        author_email = 'congcong_liu@icloud.com',
        license = 'GPLv3',
        description = 'KEGG analysis for metagenome.',
        install_requires = ['numpy',],
        scripts = ['bin/mkegg',],
        packages = find_packages(),
    )


if __name__ == '__main__':
    main()