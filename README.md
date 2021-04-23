## Motivation
When analyzing DNA sequence data of an affected individual, knowing which 
nucleotide was inherited from which parent can be beneficial when trying to 
identify certain types of DNA variants (e.g. compound heterozygous variants). 
When DNA sequence data is available for the affected individual and both 
parents (a trio), Mendelian inheritance logic can be used to determine phase. 
For example, if one parent is homozygotic for the reference allele (A/A), the 
other parent is heterozygotic (A/G), and the affected individual is 
heterozygotic (A/G), then it can be determined with surety that the reference 
allele (A) was inherited from one parent and the variant allele (G) was 
inherited from the other. However, when both parents and the affected 
individual are heterozygotic (A/G), Mendelian inheritance alone cannot be used 
to determine which nucleotide came from which parent. Phasing programs such as 
SHAPEIT4 can be used for such scenarios. SHAPEIT4 uses haplotype reference 
panels in conjunction with statistical models to help phase the data. However, 
computational phasing is not 100% accurate and phasing can not be conducted for 
positions that have no reference panel information, or positions that are 
multi-allelic. To help overcome the limitations of reference panel-based 
phasing we created trioPhaser. trioPhaser initially phases the data using 
Mendelian inheritance logic, then uses SHAPEIT4 to phase the positions that 
were not phasable using Mendelian inheritance alone. trioPhaser relies on gVCF 
files from the affected individual and his/her parents as initial input and 
outputs phased data for the affected individual. trioPhaser has been shown to
increase the total number of correctly phased positions by using Mendelian
inheritance to phase positions that SHAPEIT4 cannot phase, and to correct 
erroneously phased positions output by SHAPEIT4 [insert manuscript doi here].

## Execution
trioPhaser is executed within a Docker container, therefore, the [Docker engine
needs to be installed](https://docs.docker.com/desktop/) and the container 
needs to be downloaded 
([triophaser](https://hub.docker.com/repository/docker/dmill903/triophaser])). 
The triophaser container contains all the software that the "trio_phaser.py" 
script relies on. The "trio_phaser.py" script found in this repository is also 
contained within the container at the root directory. Therefore, when using the
triophaser Docker container, this repository (dmiller903/trioPhaser) does not 
need to be cloned. 

After installing the Docker engine and downloading the triophaser container,
"docker run -d -v {directory on local computer}:{the name you want this
directory to have inside the container} -w {the working directory within the
container} -t triophaser" is the docker code needed 
to execute the container. Everything after this bit of code is being executed 
within and by the container. The "-d" option allows the container to run in 
"detached" mode; executing the container in the background. When the "-d" 
option is used, the container ID is output upon execution. The container ID is 
output to the terminal, or so you don't have to keep track of the ID, ">" can 
be used to store the container ID to a file (as done in the example below). 
On our machine, each trio directory was stored within the "/Data" directory. 
Therefore, we attached the "/Data" directory to the container using the "-v" 
option. This allows the container to access all directories and files within 
the attached directory. We called this directory "/proj" within the container 
using ":". We then set the "/proj" directory as our working directory within 
the container using "-w". "-t" was used to allow the container to execute 
commands. "triophaser" is the name of the container to execute.
The "trio_phaser.py" script is executed within the container and has 6 required 
arguments. The arguments are required to be listed in order as follows:

1. Input file of affected individual.
2. Input file of the father.
3. Input file of the mother.
4. Name of output file.
5. Path were haplotype reference files will be saved to.
6. Number of cores available for use (up to 22, default is 2).

Here we provide example code on how to execute trioPhaser. This example implies
that the gVCF files are found within the "/Data" directory of the local machine,
that within the "/Data" directory there is a directory called
"haplotype_references" where trioPhaser will store the haplotype reference 
files that are used as part of the phasing process, and that 22 CPU cores are 
available for phasing:

```ignore
docker run -d -v /Data:/proj -w /proj -t triophaser python3 /trio_phaser.py \
  affected_individual.g.vcf.gz \
  father_of_affected_individual.g.vcf.gz \
  mother_of_affected_individual.g.vcf.gz \
  phased_output.vcf.gz \
  haplotype_references \
  22 \
  > trio_phaser.out
```