# trioPhaser: Using Mendelian inheritance logic to improve genomic phasing of trios.
Mendelian inheritance logic can be used to determine phase for the majority 
(~69%) of variant nucleotide positions. However, when both parents and the 
child are heterozygotic (A/G), Mendelian inheritance alone cannot be used to
determine phase. When such scenarios are present, phasing programs that rely on
population-based haplotype reference panels can be used.

Available software that directly uses trios and Mendelian inheritance logic 
during phasing have limitations. For example, some software require multiple 
input sources (resulting in extensive storage requirements), or do not support 
current genome builds. To address these challenges, we have developed 
trioPhaser. 

trioPhaser can increase the total number of correctly phased positions by using
Mendelian inheritance logic and *SHAPEIT4* [insert manuscript doi here].

# Usage
![trioPhaser usage](https://drive.google.com/file/d/1JM_LIaA9Mp4DeqiH7UCkTF3KfcIIHIWc/view)
## Pull the trioPhaser image (requires Docker installation):
```
docker pull dmill903/triophaser:1.0
```
## Run trioPhaser on the included example data (see [here](https://github.com/dmiller903/trioPhaser/blob/main/\validate/validate.pdf) for how we generated the files and for more details on how to run trioPhaser):

```ignore
docker run -d -v /Data:/proj -w /proj -t dmill903/triophaser:1.0 \
  python3 /trio_phaser.py \
  /trioPhaser/validate/son_GRCh38_chr22.g.vcf.gz \
  /trioPhaser/validate/father_GRCh38_chr22.g.vcf.gz \
  /trioPhaser/validate/mother_GRCh38_chr22.g.vcf.gz \
  /proj/phased_output.vcf.gz \
  /proj/haplotype_references \
  > /Data/trio_phaser.out
```
*Note: The first time trioPhaser is run, it will download reference files, which 
may take a while depending on internet speed. Subsequent runs can 
use the cached reference files.*

## View log information output by trioPhaser
```ignore
CONTAINERID="cat /Data/trio_phaser.out"
docker logs `$CONTAINERID` >> /Data/trio_phaser.out
```