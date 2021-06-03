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
parser.add_argument('haplotype_reference_files', help='The path where the \
    haplotype reference files will be/were downloaded to.\
    When using Docker, this path must be accessible by the container. If the \
    folder you want these files downloaded to are not within the same path as \
    you input files and/or output_files, you need to attach another volume to \
    the Docker container so the folder can be accessed.')
parser.add_argument('--number_of_tasks', help='Max number of cores that can be \
    used in a given run. If 22 cores are availabe for use, all chromosomes will \
    be phased by SHAPEIT4 at the same time, significantly speeding up the overall \
    runtime.', default="2")
parser.add_argument('--build_version', help='GRCh 37 or 38 are supported. GRCh38\
    is used as default.', default='38')

args = parser.parse_args()

# Create variables of each argument from argparse
child_file = args.child_file
paternal_file = args.paternal_file
maternal_file = args.maternal_file
output_file = args.output_file
haplotype_path = args.haplotype_reference_files
if not haplotype_path.endswith("/"):
    haplotype_path = haplotype_path + "/"
number_tasks = int(args.number_of_tasks)
try:
    if "38" in args.build_version:
        build_version = 38
    elif "37" in args.build_version:
        build_version = 37
except:
    print("Build version unknown. Please indicate whether input files were \
        GRCh build 37 or 38 using the --build_version argument.")

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

def get_phase(child_allele_1, child_allele_2, 
            paternal_genotype_or_haplotype, maternal_genotype_or_haplotype):
    """Function is used to get the phase of the child, using the maternal and
    paternal genotypes as a reference. The if and elif section phases the child 
    using Mendelian inheritance and lists paternal nucleotide first, and the 
    maternal nucleotide second.
    """
    if child_allele_1 in paternal_genotype_or_haplotype \
        and child_allele_1 not in maternal_genotype_or_haplotype \
        and ("." not in maternal_genotype_or_haplotype \
        and "." not in paternal_genotype_or_haplotype):
        phase = f"{child_allele_1}|{child_allele_2}"
    elif child_allele_2 in paternal_genotype_or_haplotype \
        and child_allele_2 not in maternal_genotype_or_haplotype \
        and ("." not in maternal_genotype_or_haplotype \
        and "." not in paternal_genotype_or_haplotype):
        phase = f"{child_allele_2}|{child_allele_1}"
    elif child_allele_1 in maternal_genotype_or_haplotype \
        and child_allele_1 not in paternal_genotype_or_haplotype \
        and ("." not in maternal_genotype_or_haplotype \
        and "." not in paternal_genotype_or_haplotype):
        phase = f"{child_allele_2}|{child_allele_1}"
    elif child_allele_2 in maternal_genotype_or_haplotype \
        and child_allele_2 not in paternal_genotype_or_haplotype \
        and ("." not in maternal_genotype_or_haplotype \
        and "." not in paternal_genotype_or_haplotype):
        phase = f"{child_allele_1}|{child_allele_2}"
    elif child_allele_1 == child_allele_2 \
        and (child_allele_1 in maternal_genotype_or_haplotype \
        and child_allele_1 in paternal_genotype_or_haplotype):
        phase = f"{child_allele_1}|{child_allele_2}"
    else:
        phase = "."
    return(phase)

# Create a dictionary where the key is the family members title, and the value is the sample's ID in the VCF.
sample_ids = {} #The dictionary that relate_sample_name_to_file() will use
relate_sample_name_to_file(child_file, "child")
relate_sample_name_to_file(paternal_file, "paternal")
relate_sample_name_to_file(maternal_file, "maternal")

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
temp_genotyped_name = "/tmp/genotyped.vcf.gz"

# Extract fasta reference file
os.system("unzip /fasta_references.zip -d /fasta_references")
os.system("gzip -d /fasta_references/*.gz")

# Use GATK to combine the trio into a single file and then create a gVCF
file_string = ""
for file in files:
    file_string += f"-V {file} "
    os.system(f"gatk IndexFeatureFile -F {file}")

if build_version == 38:
    os.system(f"gatk CombineGVCFs -R /fasta_references/Homo_sapiens_assembly38.fasta {file_string} -O {temp_combined_name}")
    print("Trio has been combined and written to a temporary file.")
    os.system(f"gatk IndexFeatureFile -F {temp_combined_name}")
    os.system(f"gatk --java-options '-Xmx4g' GenotypeGVCFs -R /fasta_references/Homo_sapiens_assembly38.fasta -V {temp_combined_name} -O {temp_genotyped_name}")
    print("Trio has been joint-genotyped.")
elif build_version == 37:
    os.system(f"gatk CombineGVCFs -R /fasta_references/human_g1k_v37_modified.fasta {file_string} -O {temp_combined_name}")
    print("Trio has been combined and written to a temporary file.")
    os.system(f"gatk IndexFeatureFile -F {temp_combined_name}")
    os.system(f"gatk --java-options '-Xmx4g' GenotypeGVCFs -R /fasta_references/human_g1k_v37_modified.fasta -V {temp_combined_name} -O {temp_genotyped_name}")
    print("Trio has been joint-genotyped.")

# Separate combined trio file by chromosome and create child scaffold 
# (phased VCF) for each chromosome
with gzip.open(temp_genotyped_name, "rt") as vcf:
    output_name = f"/tmp/genotyped"
    chromosome_set = set()
    header_chromosome = ""
    header_scaffold = ""
    file_list_to_bgzip = []
    for line in vcf:
        if line.startswith("##"):
            header_chromosome = header_chromosome + line
            header_scaffold = header_scaffold + line
        elif line.startswith("#CHROM"):
            header_chromosome = header_chromosome + line
            line_list = line.rstrip("\n").split("\t")
            chrom_index = line_list.index("#CHROM")
            child_index = line_list.index(sample_ids["child"])
            paternal_index = line_list.index(sample_ids["paternal"])
            maternal_index = line_list.index(sample_ids["maternal"])
            # Recreate the header line with all columns except parental columns
            new_line_list = []
            for i in range(0, len(line_list)):
                if i != maternal_index or i != paternal_index:
                    new_line_list.append(line_list[i])
            line_list = new_line_list
            line = "\t".join(line_list) + "\n"
            header_scaffold = header_scaffold + line
        elif not line.startswith("#") and line.split("\t")[0] not in chromosome_set:
            line_list = line.rstrip("\n").split("\t")
            chrom = line_list[chrom_index]
            # Chromosomes are listed as "chr1", this removes the "chr"
            updated_chrom = chrom[3:]
            line = line.replace(chrom, updated_chrom)
            # output metadata to output chromosome and output scaffold
            with gzip.open(f"{output_name}_{chrom}.vcf", "wb") as chromosome, \
                gzip.open(f"{output_name}_{chrom}_scaffold.vcf", "wb") as scaffold:
                chromosome.write(header_chromosome.encode())
                chromosome.write(line.encode())
                new_line_list = []
                for i in range(0, len(line_list)):
                    if i != maternal_index or i != paternal_index:
                        new_line_list.append(line_list[i])
                line_list = new_line_list
                line = "\t".join(line_list) + "\n"
                line = line.replace(chrom, updated_chrom)
                scaffold.write(header_scaffold.encode())
                scaffold.write(line.encode())
                chromosome_set.add(chrom)
                file_list_to_bgzip.append(f"{output_name}_{chrom}.vcf")
                file_list_to_bgzip.append(f"{output_name}_{chrom}_scaffold.vcf")
        elif not line.startswith("#") and line.split("\t")[0] in chromosome_set:
            line_list = line.rstrip("\n").split("\t")
            chrom = line_list[chrom_index]
            updated_chrom = chrom[3:]
            child_genotype = line_list[child_index].split(":")[0]
            paternal_genotype = line_list[paternal_index].split(":")[0]
            maternal_genotype = line_list[maternal_index].split(":")[0]
            child_allele_1 = child_genotype[0]
            child_allele_2 = child_genotype[-1]
            line = line.replace(chrom, updated_chrom)
            # Output line to chromosome file
            with gzip.open(f"{output_name}_{chrom}.vcf", "ab") as chromosome:
                if "." not in child_genotype and "." not in paternal_genotype \
                and "." not in maternal_genotype and (child_genotype != "0/0" \
                    or child_genotype != "0|0"):
                    chromosome.write(line.encode())
            # Phase child using Mendelian inheritance, output haplotype to 
            # scaffold file
            with gzip.open(f"{output_name}_{chrom}_scaffold.vcf", "ab") as scaffold:
                new_line_list = []
                for i in range(0, len(line_list)):
                    if i != maternal_index or i != paternal_index:
                        new_line_list.append(line_list[i])
                line_list = new_line_list
                line = "\t".join(line_list) + "\n"
                line = line.replace(chrom, updated_chrom)

                phase = get_phase(child_allele_1, child_allele_2, 
                                paternal_genotype, maternal_genotype)
                # Outputs phasable positions to the output scaffold
                if phase != ".":
                    line = line.replace(child_genotype, phase)
                    scaffold.write(line.encode())

# bgZip and index scaffold and chromosome files
for file in file_list_to_bgzip:
    bgzip_file(file)
    os.system(f"tabix -fp vcf {file}.gz")
    os.system(f"bcftools index {file}.gz")
    
print(
    "Trio has been separated into chromosome files and chromosome "
    "scaffolds have been created in preparation for phasing.")

# Extract genetic maps and download haplotype references if necessary
os.system("tar -xf /shapeit4/maps/genetic_maps.b38.tar.gz -C /shapeit4/maps/")
if build_version == 38 \
    and not os.path.exists(f"{haplotype_path}ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"):
    os.system("wget --no-check-certificate https://files.osf.io/v1/resources/rbzma/providers/osfstorage/608b8dd719183d00cb5556c3/?zip= -O /haplotype_references.zip")
    os.system(f"unzip /haplotype_references.zip -d {haplotype_path}")
    os.system(f"chmod 777 {haplotype_path}*")
    os.system("rm /haplotype_references.zip")

elif build_version == 37 \
    and not os.path.exists(f"{haplotype_path}ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"):
    os.system("wget --no-check-certificate https://files.osf.io/v1/resources/rbzma/providers/osfstorage/60b6c4639096b7023a63c8d0/?zip= -O /haplotype_references.zip")
    os.system(f"unzip /haplotype_references.zip -d {haplotype_path}")
    os.system(f"chmod 777 {haplotype_path}*")
    os.system("rm /haplotype_references.zip")

# Create a list of shapeit4 execution commands.
task_list = []
if build_version == 38:
    for i in range(22, 0, -1):
        if os.path.exists(f"/tmp/genotyped_chr{i}.vcf.gz"):
            task_list.append(f"shapeit4 --input /tmp/genotyped_chr{i}.vcf.gz --map /shapeit4/maps/chr{i}.b38.gmap.gz --region {i} --output /tmp/phased_chr{i}_with_scaffold.vcf.gz --reference {haplotype_path}ALL.chr{i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --sequencing --scaffold /tmp/genotyped_chr{i}_scaffold.vcf.gz --seed 123456789")
elif build_version == 37:
    for i in range(22, 0, -1):
        if os.path.exists(f"/tmp/genotyped_chr{i}.vcf.gz"):
            task_list.append(f"shapeit4 --input /tmp/genotyped_chr{i}.vcf.gz --map /shapeit4/maps/chr{i}.b37.gmap.gz --region {i} --output /tmp/phased_chr{i}_with_scaffold.vcf.gz --reference {haplotype_path}ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --sequencing --scaffold /tmp/genotyped_chr{i}_scaffold.vcf.gz --seed 123456789")

# Phase with shapeit4, using concurrent.futures to phase all chromosomes at once.
with concurrent.futures.ProcessPoolExecutor(max_workers=number_tasks) as executor:
    executor.map(os_system_task, task_list)

# Iterate through the phased file and determine if SHAPEIT4 phased the 
# positions correctly that can be phased using Mendelian inheritance.
shapeit_positions = {}
correctly_phased = 0
incorrectly_phased = 0
total_variants = 0
total_phased = 0
could_not_be_determined = 0
header = ""
header_written = False
with gzip.open(output_file.replace('.vcf.gz', \
    '_phased_incorrectly_by_SHAPEIT4.vcf.gz'), "wb") as incorrectly_phased_out:
    for i in range(1, 23):
        if os.path.exists(f"/tmp/phased_chr{i}_with_scaffold.vcf.gz"):
            with gzip.open(f"/tmp/phased_chr{i}_with_scaffold.vcf.gz", "rt") as phasedFile:
                for line in phasedFile:
                    if line.startswith("##"):
                        if i == 1:
                            header = header + line
                    elif line.startswith("#CHROM"):
                        line_list = line.rstrip("\n").split("\t")
                        child_index = line_list.index(sample_ids["child"])
                        paternal_index = line_list.index(sample_ids["paternal"])
                        maternal_index = line_list.index(sample_ids["maternal"])
                        pos_index = line_list.index("POS")
                        chrom_index = line_list.index("#CHROM")
                        info_index = line_list.index("INFO")
                        if header_written is False:
                            header = header + line
                            incorrectly_phased_out.write(header.encode())
                            header_written = True
                    else:
                        total_phased += 1
                        total_variants += 1
                        line_list = line.rstrip("\n").split("\t")
                        chrom = line_list[chrom_index]
                        pos = int(line_list[pos_index])
                        child_haplotype = line_list[child_index]
                        paternal_haplotype = line_list[paternal_index]
                        maternal_haplotype = line_list[maternal_index]
                        child_allele_1 = child_haplotype[0]
                        child_allele_2 = child_haplotype[-1]
                        line_list[info_index] = "."
                        line_list[chrom_index] = "chr" + chrom

                        if chrom not in shapeit_positions:
                            shapeit_positions[chrom] = {}
                        phase = get_phase(child_allele_1, child_allele_2, 
                                        paternal_haplotype, maternal_haplotype)

                        # Checks if the "child_haplotype" (phased with SHAPEIT4) was
                        # correctly phased using Mendelian determined "phase".
                        if phase == child_haplotype and phase != ".":
                            correctly_phased += 1
                            line = "\t".join(line_list) + "\n"
                            shapeit_positions[chrom][pos] = line
                        elif phase != child_haplotype and phase != ".":
                            incorrectly_phased += 1
                            incorrectly_phased_out.write(line.encode())
                            line_list[child_index] = phase
                            line = "\t".join(line_list) + "\n"
                            shapeit_positions[chrom][pos] = line
                        else:
                            could_not_be_determined += 1
                            line = "\t".join(line_list) + "\n"
                            shapeit_positions[chrom][pos] = line

# Iterate through the genotyped file and if the position was not phased by 
# SHAPEIT4, determine if it's phasable and if it is, put it in the final output.
not_in_shapeit = 0
for i in range(1, 23):
    if os.path.exists(f"/tmp/genotyped_chr{i}.vcf.gz"):
        with gzip.open(f"/tmp/genotyped_chr{i}.vcf.gz", "rt") as genotypeFile:
            for line in genotypeFile:
                if line.startswith("##"):
                    continue
                elif line.startswith("#CHROM"):
                    line_list = line.rstrip("\n").split("\t")
                    child_index = line_list.index(sample_ids["child"])
                    paternal_index = line_list.index(sample_ids["paternal"])
                    maternal_index = line_list.index(sample_ids["maternal"])
                    pos_index = line_list.index("POS")
                    chrom_index = line_list.index("#CHROM")
                    filter_index = line_list.index("FILTER")
                    info_index = line_list.index("INFO")
                    format_index = line_list.index("FORMAT")
                else:
                    line_list = line.rstrip("\n").split("\t")
                    chrom = line_list[chrom_index]
                    pos = int(line_list[pos_index])
                    child_haplotype = line_list[child_index].split(":")[0]
                    paternal_genotype = line_list[paternal_index].split(":")[0]
                    maternal_genotype = line_list[maternal_index].split(":")[0]
                    child_allele_1 = child_haplotype[0]
                    child_allele_2 = child_haplotype[-1]
                    # If the position is not in the shapeit_positions dictionary, 
                    # then Mendelian inheritance logic is used to see if position 
                    # is phasable.
                    if pos not in shapeit_positions[chrom]:
                        total_variants += 1
                        line_list[filter_index] = "."
                        line_list[info_index] = "."
                        line_list[format_index] = "GT"
                        line_list[chrom_index] = "chr" + chrom
                        line_list[paternal_index] = paternal_genotype
                        line_list[maternal_index] = maternal_genotype
                        # Phase using Mendelian inheritance when able.
                        phase = get_phase(child_allele_1, child_allele_2, 
                                        paternal_genotype, maternal_genotype)
                        if phase != ".":
                            not_in_shapeit += 1
                            total_phased += 1
                            line_list[child_index] = phase
                            line = "\t".join(line_list) + "\n"
                            shapeit_positions[chrom][pos] = line

# Print summary statistics.
print(f"\nThere were {correctly_phased} ({(correctly_phased / (correctly_phased + incorrectly_phased)) * 100:.5f}%) correctly phased haplotypes, and {incorrectly_phased} ({(incorrectly_phased / (correctly_phased + incorrectly_phased)) * 100:.5f}%) incorrectly phased haplotypes. (as phased by SHAPEIT4)")
print(f"The {incorrectly_phased} ({(incorrectly_phased / (correctly_phased + incorrectly_phased)) * 100:.5f}%) incorrectly phased haplotypes were corrected using Mendelian inheritance and will be included the phased output file")
print(f"{not_in_shapeit} variants were not phased by shapeit but were phasable and will be included in final output.")
print(f'{could_not_be_determined} phased variants (as phased by SHAPEIT4) could not be verified by using Mendelian inheritance alone.')
print(f"There were {total_phased} total variants phased.\n")
print(f"There were {total_variants} total variants prior to any phasing.\n")

# Output the final file.
with gzip.open(output_file.replace(".gz", ""), "wb") as output:
    output.write(header.encode())
    for chrom, posDict in sorted(shapeit_positions.items()):
        for pos, line in sorted(posDict.items()):
            output.write(line.encode())

# bgZip and index final file.
os.system(f"zcat {output_file.replace('.gz', '')} | bgzip -f > {output_file}")
os.system(f"tabix -fp vcf {output_file}")
os.system(f"bcftools index {output_file}")
os.system(f"chmod 777 {output_file} {output_file}.tbi {output_file}.csi")
os.system(f"rm {output_file.replace('.gz', '')}")
print(f"\nPhased output file written as {output_file}\n")

# Print message and how long the previous steps took.
timeElapsedMinutes = round((time.time()-start_time) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours) {char}')