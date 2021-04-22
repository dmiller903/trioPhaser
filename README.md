# trioPhaser
When analyzing DNA sequence data of an affected individual, knowing which 
nucleotide was inherited from which parent can be beneficial when trying to 
uncover certain types of DNA variants (e.g. compound heterozygous variants). 
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