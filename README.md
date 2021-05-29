# RNPFind

## What
RNPFind is a bioinformatics tool. It takes a gene name as input and gives a
list of RNA binding proteins (RBPs) that bind to the input gene's RNA
transcript and their binding sites on the transcript as an output. In addition
to collecting the data needed to produce this output, RNPFind provides options
for analysing the data collected (for example, suggesting associations among
the RBPs).

## When
This tool was produced at iYLab at CMUQ to help both my molecular biology
research project as well as that of others in the lab. Originally a CLI tool
that needed to be installed on the user's computer, it was made into a webapp
later for ease of use and access.

## Where
The webapp is avaiable for use at https://rnp-find.herokuapp.com

## How
RNPFind collects its data from three main databases: ATTRACT, RBPDB, and
POSTAR. Each database stores binding sites of RBPs deduced through various
methods, such as various experimental methods as well as computational methods.
RBPDB and ATTRACT provide binding site patterns for RBPs, so we scan the input
gene's sequence to deduce the binding sites. POSTAR stores experimentally
deduced binding sites for the transcriptome for all RBPs, so we simply query.

The various analysis methods are detailed elsewhere.

## Steps for Running
For now, `rnpfind` requires my AWS credentials to fetch data. This will be
fixed soon, but for now the following steps are just for me.

Load AWS credentials into the environment
```bash
source env/bin/activate
```

*Note:* In the following commands, docker might need `sudo` prepended depending
on your setup.

Build the Docker image if needed
```bash
docker build --tag rnpfind-test .
```

Run it on any port (here we run on port 80)
```bash
docker run -p 80:8000 -e PORT=8000 -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
-e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY --name rnpfind rnpfind-test
```

or just run `./start.sh` (which does the same)

Now you should be able to use `rnpfind` by navigating to localhost:80 (or the
public IP of your machine)

Once done, stop and remove by:
```bash
docker stop rnpfind
docker rm rnpfind
```

or `./stop.sh` (which does the same)




