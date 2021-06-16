# trioPhaser: Using Mendelian inheritance logic to improve genomic phasing of trios.
Mendelian inheritance logic can be used to determine phase for the majority 
(67-83%) of variant nucleotide positions (Choi et al. 2018 & Roach et al. 2011). 
However, when both parents and the  child are heterozygotic (A/G), Mendelian 
inheritance alone cannot be used to determine phase. When such scenarios are 
present, phasing programs that rely on population-based haplotype reference 
panels can be used.

Available software that directly uses trios and Mendelian inheritance logic 
during phasing have limitations. For example, some software require multiple 
input sources (resulting in extensive storage requirements), or do not support 
current genome builds. To address these challenges, we have developed 
trioPhaser. 

trioPhaser can increase the total number of correctly phased positions by using
Mendelian inheritance logic in conjunction with *SHAPEIT4* [insert manuscript doi here].

# Usage
<a href="https://drive.google.com/uc?export=view&id=1LXPqaC5SG05zVG8R_sCyaieGJIHeAWAC">
    <img src="https://drive.google.com/uc?export=view&id=1LXPqaC5SG05zVG8R_sCyaieGJIHeAWAC"
    style="width: 300px; max-width: 100%; height: auto"
    title="Click for the larger version." />
</a>

## Installation (requires [Docker installation](https://docs.docker.com/desktop/))
```
docker pull dmill903/triophaser:latest
```
## Run trioPhaser on the included example data (see [here](https://github.com/dmiller903/trioPhaser/blob/main/validate/validate.pdf) for how we generated the files and for more details on how to run trioPhaser)
Here we provide example code on how to execute trioPhaser. In this example, 
the input files are provided within the container. Under normal circumstances, 
the input files will be on your local machine. This example implies the following:

1. There is a directory called `/Data` on the local machine.
2. We attached and called this directory `/proj` within the container.
3. We set `/proj` as our working directory within the container.
4. The `haplotype_references` directory has been created and is located at 
`/Data/haplotype_references` (*this is where trioPhaser will store the 
haplotype reference data).
5. The output file will be written to `/Data/phased_output.vcf.gz`
6. The container ID will be stored at `/Data/trio_phaser.out` upon execution. 

Update `/Data` and `haplotype_references` with the directories you want to use 
on your local machine. For example, if you want the output to be stored at 
`/tmp` and within `/tmp` you have a folder called `references`, you would 
change `/Data` to `/tmp` and `haplotype_references` to `references`.

Use `docker run -t dmill903/triophaser:latest python3 /trio_phaser.py -h`
learn more about trioPhaser's arguments.

```ignore
docker run -d -v /Data:/proj -w /proj -t dmill903/triophaser:latest \
  python3 /trio_phaser.py \
  /trioPhaser/validate/son_GRCh38_chr22.g.vcf.gz \
  /trioPhaser/validate/father_GRCh38_chr22.g.vcf.gz \
  /trioPhaser/validate/mother_GRCh38_chr22.g.vcf.gz \
  /proj/phased_output.vcf.gz \
  /proj/haplotype_references \
  > /Data/trio_phaser.out
```

*\*Note: The first time trioPhaser is run, it will download reference files, 
which may take a while depending on internet speed. Subsequent runs can use the
cached reference files.*

## View log information output by trioPhaser
```ignore
CONTAINERID="cat /Data/trio_phaser.out"
docker logs `$CONTAINERID` >> /Data/trio_phaser.out
```

## References
Choi Y, Chan AP, Kirkness E, Telenti A, Schork NJ: Comparison of phasing 
  strategies for whole human genomes. PLoS Genet. 2018, 14:e1007308.

Roach JC, Glusman G, Hubley R, Montsaroff SZ, Holloway AK, Mauldin DE, Srivastava D, 
  Garg V, Pollard KS, Galas DJ, Hood L, Smit AFA: Chromosomal haplotypes by 
  genetic phasing of human families. Am. J. Hum. Genet. 2011, 89:382–397.