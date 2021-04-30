## Motivation
When analyzing DNA sequence data of an individual, knowing which nucleotide was 
inherited from which parent can be beneficial when trying to identify certain 
types of DNA variants (e.g. compound heterozygous variants). When DNA sequence 
data is available for the individual and both parents (a trio), Mendelian 
inheritance logic can be used to determine phase. For example, if one the 
mother is homozygotic for the reference allele (A/A), the father is 
heterozygotic (A/G), and the child is heterozygotic (A/G), then it can be 
determined that the reference allele (A) was inherited from the mom and the 
variant allele (G) was inherited from the dad. However, when both parents and 
the affected individual are heterozygotic (A/G), Mendelian inheritance alone 
cannot be used to determine which nucleotide came from which parent. Phasing 
programs such as *SHAPEIT4* can be used for such scenarios. *SHAPEIT4* uses a
haplotype reference panel in conjunction with statistical models to help phase 
the data. However, computational phasing is not 100% accurate and phasing can 
not be conducted for positions that have no reference panel information, or 
positions that are multi-allelic. To help overcome the limitations of reference
panel-based phasing we created trioPhaser. trioPhaser initially phases the data
using Mendelian inheritance logic, then uses *SHAPEIT4* to phase the positions 
that were not phaseable using Mendelian inheritance alone. trioPhaser relies on 
gVCF files from a child and his/her parents as initial input and 
outputs a phased VCF file. trioPhaser has been shown to increase the total 
number of correctly phased positions by using Mendelian inheritance to phase 
positions that *SHAPEIT4* cannot phase, and to correct erroneously phased 
positions output by *SHAPEIT4* [insert manuscript doi here].

## Execution
trioPhaser is executed within a Docker container, therefore, the [Docker engine
needs to be installed](https://docs.docker.com/desktop/) and the container 
needs to be downloaded ([triophaser](https://hub.docker.com/repository/docker/dmill903/triophaser])). The triophaser container contains all the software that the "trio_phaser.py" script relies on.
The "trio_phaser.py" script found in this repository is also contained within 
the container at the root directory. Therefore, when using the triophaser 
Docker container, this repository (dmiller903/trioPhaser) does not need to be 
cloned. 

To download trioPhaser after docker is installed, use:
```
docker pull dmill903/triophaser:1.0
```

"docker run -d -v {directory on local computer}:{the name you want this
directory to have inside the container} -w {the working directory within the
container} -t dmill903/triophaser:1.0" is the docker code needed 
to execute the container. Everything after this bit of code is being executed 
within and by the container. The "-d" option allows the container to run in 
"detached" mode; executing the container in the background. When the "-d" 
option is used, the container ID is output upon execution. The container ID is 
output to the terminal, or so you don't have to keep track of the ID, ">" can 
be used to store the container ID to a file (as done in the example below). 
In order for the container to be able to find files on a local machine, a 
volume needs to be attached with the Docker argument "-v". The container can 
see any directories and files contained within the attached volume. For 
example, if you have a directory called "/Data" on your local machine and you
attach this directory to the container, then the container can access anything
within the "/Data" directory. After specifying the directory you want to attach
":" can be used to specify what this directory will be called within the
container. For example, if you want the "/Data" directory to be called "/proj"
within the container, you would used "-v /Data:/proj". You can specify what
directory the container will use as its working directory with "-w". For 
example, we can set the "/proj" directory as the working directory with
"-w /proj". "-t" is used to allow the container to execute commands. 
"dmill903/triophaser:1.0" is the name of the container to execute. The 
"trio_phaser.py" script is executed within the container (its located at the
root directory within the container) and has 6 required arguments. The 
arguments need to be listed in order as follows:

1. Input file of affected individual (path accessible by the attached volume).
2. Input file of the father (path accessible by the attached volume).
3. Input file of the mother (path accessible by the attached volume).
4. Name of output file (path accessible by the attached volume).
5. Path were haplotype reference files will be saved to (requires ~13GB of
storage, and needs to be accessible by the attached volume).
6. Number of cores available for use (up to 22, default is 2).

Here we provide example code on how to execute trioPhaser. In this example,
the input files are provided within the container. Under normal circumstances,
the input files will be on your local machine. Details on how the provided 
files were created can be found within our  
[validation documentation](https://github.com/dmiller903/trioPhaser/blob/main/\validate/validate.pdf). This example implies that 1) there is a directory called "/Data" on the local machine, 2) we
attached and called this directory "/proj" within the container, 3) we set 
"/proj" as our working directory within the container, 4) the 
"haplotype_references" directory has been created and is located at 
"/Data/haplotype_references" (*this is where trioPhaser will store the haplotype
reference data), 5) the output file will be written to 
"/Data/phased_output.vcf.gz", 6) the container ID will be stored at 
"/Data/trio_phaser.out" upon execution, and 7) 2 CPU cores are available on the
local machine for trioPhaser to use.

*Note: The first time trioPhaser is executed, the haplotype reference files
are downloaded to the path specified. This download process can take awhile
depending on internet speeds. Therefore, expect a longer than usual run-time. 
Upon subsequent trioPhaser runs, the reference files will not be downloaded as 
long as trioPhaser has access to the path where these files are stored.

```ignore
docker run -d -v /Data:/proj -w /proj -t dmill903/triophaser:1.0 \
  python3 /trio_phaser.py \
  /trioPhaser/validate/son_GRCh38_chr22.g.vcf.gz \
  /trioPhaser/validate/father_GRCh38_chr22.g.vcf.gz \
  /trioPhaser/validate/mother_GRCh38_chr22.g.vcf.gz \
  phased_output.vcf.gz \
  haplotype_references \
  2 \
  > /Data/trio_phaser.out
```

To look at the log file from the above example, on your local machine, type:

```ignore
CONTAINERID="cat /Data/trio_phaser.out"
docker logs `$CONTAINERID` >> /Data/trio_phaser.out
```