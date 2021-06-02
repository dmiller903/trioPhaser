# Import necessary modules
import gzip
import os
import time
import concurrent.futures
import argparse

# Keep track of when the script began
start_time = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Phases a trio when gVCF files are \
    available for each individual. This program requires 14 GB of haplotype \
    reference files and these files will automatically be downloaded when the \
    program is first executed.')
parser.add_argument('child_file', help='Sample (patient) File. Must be gzipped')
parser.add_argument('paternal_file', help='Paternal File. Must be gzipped')
parser.add_argument('maternal_file', help='Maternal File. Must be gzipped')
parser.add_argument('output_file', help='Name and path of output file')
parser.add_argument('--number_of_tasks', help='Max number of cores that can be \
    used in a given run. If 22 cores are availabe for use, all chromosomes will \
    be phased by SHAPEIT4 at the same time, significantly speeding up the overall \
    runtime.', default="2")

args = parser.parse_args()

# Create variables of each argument from argparse
child_file = args.child_file
paternal_file = args.paternal_file
maternal_file = args.maternal_file
output_file = args.output_file
number_tasks = int(args.number_of_tasks)

# Functions
def relate_sample_name_to_file(file, title):
    """This function gives each input file a title of "child", "maternal", or
    paternal.
    """
    with gzip.open(file, 'rt') as gVCF:
        for line in gVCF:
            if line.startswith('##'):
                continue
            elif line.startswith("#CHROM"):
                line_list = line.rstrip("\n").split("\t")
                sample_id = line_list[-1]
                sample_ids[title] = sample_id
                break

def bgzip_file(file):
    os.system(f"zcat {file} | bgzip -f > {file}.gz")
    os.system(f"rm {file}")

def filter_child(file):
    """Filter child file, remove  variants-only sites, and create a dictionary 
    of variant-only sites.
    """
    temp_output = "/tmp/child_parsed.vcf"
    with gzip.open(file, 'rt') as gVCF, gzip.open(temp_output, 'wb') as parsed:
        for line in gVCF:
            if line.startswith('##'):
                parsed.write(line.encode())
            elif line.startswith("#CHROM"):
                line_list = line.rstrip("\n").split("\t")
                chrom_index = line_list.index("#CHROM")
                pos_index = line_list.index("POS")
                parsed.write(line.encode())
            elif "END" not in line:
                line_list = line.rstrip("\n").split("\t")
                chrom = line_list[chrom_index]
                pos = line_list[pos_index]
                if chrom not in position_dict and chrom[3:].isnumeric() and int(chrom[3:]) in range(1, 23):
                    position_dict[chrom] = {pos}
                    parsed.write(line.encode())
                elif chrom in position_dict:
                    position_dict[chrom].add(pos)
                    parsed.write(line.encode())
    bgzip_file(temp_output)

def filter_parents(file):
    """Filter each parent file for sites that occur as variants in the child of 
    that family.
    """
    if file == paternal_file:
        temp_output = f"/tmp/paternal_parsed.vcf"
    elif file == maternal_file:
        temp_output = f"/tmp/maternal_parsed.vcf"
    with gzip.open(file, 'rt') as gVCF, gzip.open(temp_output, 'wb') as parsed:
        for line in gVCF:
            if line.startswith("#"):
                parsed.write(line.encode())
            else:
                line_list = line.split("\t")
                chrom = line_list[0]
                pos = line_list[1]
                if chrom in position_dict and pos in position_dict[chrom]:
                    parsed.write(line.encode())
                elif chrom in position_dict and pos not in position_dict[chrom]:
                    if "END" in line:
                        for i in range(int(pos), int(line_list[7].lstrip("END=")) + 1):
                            if str(i) in position_dict[chrom]:
                                parsed.write(line.encode())
    bgzip_file(temp_output)
    print(f"Positions in {file} that correspond to variant-only positions of child have been output to temporary file.")

def os_system_task(task):
    """Function is used to download files or phase at various times throughout
    the process.
    """
    os.system(task)

# Create a dictionary where the key is the family members title, and the value is the sample's ID in the VCF.
sample_ids = {} #The dictionary that relate_sample_name_to_file() will use
relate_sample_name_to_file(child_file, "child")
relate_sample_name_to_file(paternal_file, "paternal")
relate_sample_name_to_file(maternal_file, "maternal")
print(sample_ids)

# Create a dictionary that has all the variant positions of the child, for each 
# chromosome and output a file that has variant-only positions for the child
position_dict = {} #The dictionary that filter_child() will use
filter_child(child_file)
print("Variant-only positions of the child have been written to a temporary file.")

# Output a temporary file for each parent that has positions that occur as variant-only positions in the child
with concurrent.futures.ProcessPoolExecutor(max_workers=number_tasks) as executor:
    executor.map(filter_parents, [paternal_file, maternal_file])

# Use GATK to combine all trios into one temporary vcf and then genotype the combined trio vcf
files = ["/tmp/child_parsed.vcf.gz", "/tmp/paternal_parsed.vcf.gz", "/tmp/maternal_parsed.vcf.gz"]
temp_combined_name = "/tmp/combined.vcf.gz"
try:
    file_string = ""
    for file in files:
        file_string += f"-V {file} "
        os.system(f"gatk IndexFeatureFile -F {file}")
    # Extract fasta reference file
    os.system("unzip /fasta_references.zip -d /fasta_references")
    os.system("gzip -d /fasta_references/*.gz")
    os.system(f"gatk CombineGVCFs -R /fasta_references/Homo_sapiens_assembly38.fasta {file_string} -O {temp_combined_name}")
    print("Trio has been combined and written to a temporary file.")
    os.system(f"gatk IndexFeatureFile -F {temp_combined_name}")
    os.system(f"gatk --java-options '-Xmx4g' GenotypeGVCFs -R /fasta_references/Homo_sapiens_assembly38.fasta -V {temp_combined_name} -O {output_file}")
    print("Trio has been joint-genotyped.")

except:
    print("Trio not combined, there was an error detected by GATK")