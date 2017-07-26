[![PyPI](https://img.shields.io/pypi/v/simulatesv.svg)](https://pypi.python.org/pypi/simulatesv/)
[![PyPI](https://img.shields.io/pypi/l/simulatesv.svg)](https://pypi.python.org/pypi/simulatesv/)

simulatesv
===

Simulate Structural Variations and SNPs for artificial DNA templates. The DNA
templates are generated simply with a pseudo-random number generator, and is
not guided by biological significance. The tool is provided with the intention
of creating a validation set of variants to test current structural variation
detection tools.

## Installation

The recommended method is to install via pip or git clone. You can also get the
source distribution from [pypi](https://pypi.python.org/pypi/simulatesv/).

`$ pip install simulatesv`

or

`$ git clone git@github.com:mlliou112/simulatesv.git`


## Usage

This script runs with either Python 2 or 3, and requires the Python package `numpy` to run.

If you install via pip, the script will be available to run as `simulatesv.py`, or `python simulatesv.py` if you installed via git clone.

#### Basic

```
$ python simulatesv.py
```
This will run the program with all defaults, making one 50 kbp reference genome and three simulated genomes in the current working directory. (total 7 files)

```
.
├── reference.fna  # randomly generated dna template on which mutations are made
├── sim_0.fna      # first simulated genome based off reference with SNPs and SVs
├── sim_1.fna      # 2nd ...
├── sim_2.fna      # 3rd ...
├── changes_0.txt  # log of what changes were made to the reference genome to create the first simulated genome
├── changes_1.txt  # 2nd change log ...
└── changes_2.txt  # 3rd change log ...

```

If you have your own genome template that you would like to simulate variants
from, you can specify the FASTA file with the `-t FILE` option. You can also
modify the frequency of each variant as well as their size. See the Command
Reference below for more options.

```
$ python simulatesv.py -t reference.fna --snp-error-rate .005 --largest-insertion-size 200
```

#### Changes


```
type       ref_idx	alt_idx	size	ref	alt
SNP	    22	     22	     1	   A	  G
INS	    33272	  .	      17	  .      CAACTACTAATCCACCA
```

The changes file logs what changes are made to the reference sequence. A `"."`
marks that column as not applicable.

| column  | description |
|---------|:----------------------------------------|
| type    | type of structural variation           |
| ref_idx | index in reference of ref bases        |
| alt_idx | index in simulated genome of alt bases |
| ref     | bases in reference genome              |
| alt     | bases in simulated genome              |

#### Notes

* If you are simulating a small DNA template (< 500), make sure that SV option
sizes are less than the total size of your genomes, otherwise you may run into
unknown issues.

* Mutation rates for SVs and SNPs are approximate. A binomial model is used to determine the actual number of mutations based off of the given error rate.

## Command Reference 

```
simulate.py [-h] [-n N] [-l N | -t FILE] [-b BASE] [-o FILE] [-c, BASE]
            [-se ERR] [-ie ERR] [-de ERR] [-lis N] [-lds N] [-sis N]
            [-sds N] [-te ERR] [-lts N] [-sts N] [-ce ERR] [-lcn N]
            [-scn N] [-lcs N] [-scs N]

Simulate DNA template sequences with known SNPs and structural variations.

optional arguments:
  -h, --help            show this help message and exit
  -n N, --number N      Number of simulated genome sequences to generate.
                        (default=3)
  -l N, --length N      Total length of reference genome sequence.
                        (default=50000)
  -t FILE, --template FILE
                        Template to use as a base for generating SNP/SV in
                        Fasta format.

Output File Options:
  -b BASE, --basename BASE
                        Basename of the simulated fake sequences.
                        (default=sim_)
  -o FILE, --output FILE
                        Reference file for the mutated template
                        (default=reference.fna)
  -c, BASE, --changes BASE
                        Basename of the changes files (default=changes_)

SNP Options:
  -se ERR, --snp-error-rate ERR
                        Error rate for SNPs (default=.001)

Indel Options:
  -ie ERR, --insertion-error-rate ERR
                        Error rate for insertions (default=.0001)
  -de ERR, --deletion-error-rate ERR
                        Error rate for deletions (default=.0001)
  -lis N, --largest-insertion-size N
                        Largest insertion size (default=500)
  -lds N, --largest-deletion-size N
                        Largest deletion size (default=500)
  -sis N, --smallest-insertion-size N
                        Smallest insertion size (default=1)
  -sds N, --smallest-deletion-size N
                        Smallest deletion size(default=1)

Trans Options:
  -te ERR, --trans-error-rate ERR
                        Error rate for translocations (default=.0001)
  -lts N, --largest-trans-size N
                        Largest translocation size (default=500)
  -sts N, --smallest-trans-size N
                        Smallest translocation size (default=10)

CNV Options:
  -ce ERR, --cnv-error-rate ERR
                        Error rate for CNVs (default=.0001)
  -lcn N, --largest-cnv-number N
                        Maximum number of times the CNV can repeat
  -scn N, --smallest-cnv-number N
                        Minimum number of times the CNV will repeat
  -lcs N, --largest-cnv-size N
                        Largest size for CNVs (default=300)
  -scs N, --smallest-cnv-size N
                        Smallest size for CNVs (default=2)

```


## Contributing

Bug reports, feature requests, and pull requests are encouraged!

## License
MIT
