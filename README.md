This is the main repository of the [NOMAD](https://www.nomad-coe.eu/) parser for
[VASP](http://www.vasp.at/).

# Example
```python
    from vaspparser import VASPRunParserInterface
    import matplotlib.pyplot as mpl

    # 1. Initialize a parser with a set of default units.
    default_units = ["eV"]
    parser = VASPRunParserInterface(default_units=default_units)

    # 2. Parse a file
    path = "path/to/main.file"
    results = parser.parse(path)

    # 3. Query the results with using the id's created specifically for NOMAD.
    scf_energies = results["energy_total_scf_iteration"]
    mpl.plot(scf_energies)
    mpl.show()
```

# Installation
The code is python 3 compatible. First download and install
the nomadcore package:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/python-common.git
cd python-common
pip install -r requirements.txt
pip install -e .
```

Then download the metainfo definitions to the same folder where the
'python-common' repository was cloned:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info.git
```

Finally download and install the parser:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/parser-vasp.git
cd parser-vasp
pip install -e .
```
