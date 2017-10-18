from setuptools import setup, find_packages


def main():
    setup(
        name="vaspparser",
        version="0.1",
        description="NOMAD parser implementation for VASP.",
        author="Fawzi Mohamed",
        license="",
        package_dir={'': 'parser/parser-vasp'},
        packages=find_packages(),
        install_requires=[
            'nomadcore'
        ],
    )

if __name__ == "__main__":
    main()
