# RNPFind

![rnpfind](images/rnpfind-logo-transparent.png "RNPFind - Explore RNA-Protein
interactions")

Explore RNA-Protein interactions.

RNPFind is a bioinformatics tool. It takes an RNA transcript as input and gives
a list of RNA binding protein (RBP) binding sites on the transcript as output.

Various output formats for representing the binding sites is supported.


## Usage
Use `rnpfind` as either a command line tool or through the webtool.

#### Command line tool
The full docs are available [here](cli/README.md), but a quick example for
installation and use is:

```bash
pip install rnpfind
rnpfind -f bed -o ./malat1-sites malat1
```

which installs `rnpfind` from the
[PyPi respository](https://pypi.org/project/rnpfind/), and then invokes the
program to generate binding sites for the gene `malat1` in the `bed` file format
in the directory `./malat1-sites`.


Note that the tool will download large data files upon first invocation (around
6.4 GB). If the memory footprint is too much for you, the web tool is probably
better for you.


#### Web tool
Visit the [RNPFind website](https://rnpfind.com) and the rest is pretty
self-explanatory.


## How does it work?
RNPFind collects its data from three main databases:
 - [ATTRACT](https://attract.cnic.es)
 - [RBPDB](http://rbpdb.ccbr.utoronto.ca)
 - [POSTAR](http://postar.ncrnalab.org)

Each database stores binding sites of RBPs deduced through various
methods, such as various experimental methods as well as computational methods.
RBPDB and ATTRACT provide binding site patterns for RBPs, so we scan the input
gene's sequence to deduce the binding sites. POSTAR stores experimentally
deduced binding sites for the transcriptome for all RBPs, so we simply query.




## Development
Enable recommended Git Hooks as follows:
```
git config --local core.hooksPath .githooks/
```
The above will run the following to ensure code consistency every time you
commit:
 - [black](https://github.com/psf/black)
 - [isort](https://github.com/PyCQA/isort)

Also use [fit-commit](https://github.com/m1foley/fit-commit) to ensure
consistent commit message style.


### Steps for running the website
Use Python 3.8 or higher

Create a file called `db.env` (choose a db, username, and password):
```bash
# Set PostgreSQL credential env variables
POSTGRES_DB=******
POSTGRES_USER=*******
POSTGRES_PASSWORD=*****************
```

*Note:* In the following commands, docker might need `sudo` prepended depending
on your setup.

Build the Docker image(s) if needed
```bash
docker-compose build
```

Run it
```bash
docker-compose up
```

Now you should be able to use `rnpfind` by navigating to localhost:80 (or the
public IP of your machine)

Once done, stop and remove by typing Ctrl+C or using `docker-compose down`


## Acknowledgements
This tool was produced at iYLab at CMUQ to help both my molecular biology
research project as well as that of others in the lab. Originally a CLI tool
that needed to be installed on the user's computer, it was made into a webapp
later for ease of use and access.
