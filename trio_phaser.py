# Import necessary modules
import gzip
import re
import os
import time
import concurrent.futures
import argparse
import glob

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Phases a trio when gVCF files are available for each individual. This program \
    requires 14 GB of haplotype reference files and these files will automatically be downloaded when the program is first \
    executed.')

parser.add_argument('child_file', help='Sample (patient) File. Must be gzipped')
parser.add_argument('paternal_file', help='Paternal File. Must be gzipped')
parser.add_argument('maternal_file', help='Maternal File. Must be gzipped')
parser.add_argument('output_file', help='Name and path of output file')
parser.add_argument('haplotype_reference_files', help='The path where the haplotype reference files will be downloaded to.\
    When using Docker, this path must be accessible by the container. If the folder you want these files downloaded to are not\
    within the same path as you input files and/or outputfiles, you need to attach another volume to the Docker container\
    so the folder can be accessed.')
parser.add_argument('number_of_tasks', help='max number of cores that can be used in a given run.')

args = parser.parse_args()

#Create variables of each argument from argparse
childFile = args.child_file
paternalFile = args.paternal_file
maternalFile = args.maternal_file
outputFile = args.output_file
haplotypePath = args.haplotype_reference_files
if not haplotypePath.endswith("/"):
    haplotypePath = haplotypePath + "/"
numberTasks = int(args.number_of_tasks)

#Functions
def relate_sample_name_to_file(file, title):
    with gzip.open(file, 'rt') as gVCF:
        for line in gVCF:
            if line.startswith('##'):
                continue
            elif line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                sampleId = lineList[-1]
                sampleIds[title] = sampleId
                break

def bgzipFile(file):
    os.system(f"zcat {file} | bgzip -f > {file}.gz")
    os.system(f"rm {file}")

#Filter child file, remove  variants-only sites, create a dictionary of variant-only sites
def filterChild(file):
    tempOutput = "/tmp/child_parsed.vcf"
    with gzip.open(file, 'rt') as gVCF, gzip.open(tempOutput, 'wb') as parsed:
        for line in gVCF:
            if line.startswith('##'):
                parsed.write(line.encode())
            elif line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                chromIndex = lineList.index("#CHROM")
                posIndex = lineList.index("POS")
                parsed.write(line.encode())
            elif "END" not in line:
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex]
                pos = lineList[posIndex]
                if chrom not in positionDict and chrom[3:].isnumeric() and int(chrom[3:]) in range(1, 23):
                    positionDict[chrom] = {pos}
                    parsed.write(line.encode())
                elif chrom in positionDict:
                    positionDict[chrom].add(pos)
                    parsed.write(line.encode())
    bgzipFile(tempOutput)

#Filter each parent file for sites that occur in sample of that family
def filterParents(file):
    if file == paternalFile:
        tempOutput = f"/tmp/paternal_parsed.vcf"
    elif file == maternalFile:
        tempOutput = f"/tmp/maternal_parsed.vcf"
    with gzip.open(file, 'rt') as gVCF, gzip.open(tempOutput, 'wb') as parsed:
        for line in gVCF:
            if line.startswith("#"):
                parsed.write(line.encode())
            else:
                lineList = line.split("\t")
                chrom = lineList[0]
                pos = lineList[1]
                if chrom in positionDict and pos in positionDict[chrom]:
                    parsed.write(line.encode())
                elif chrom in positionDict and pos not in positionDict[chrom]:
                    if "END" in line:
                        for i in range(int(pos), int(lineList[7].lstrip("END=")) + 1):
                            if str(i) in positionDict[chrom]:
                                parsed.write(line.encode())
    bgzipFile(tempOutput)
    print(f"Positions in {file} that correspond to variant-only positions of child have been output to temporary file.")

#Create a function to use to download files or phase
def osSystemTask(task):
    os.system(task)

#Create a dictionary where the key is the family members title, and the value is the sample's ID in the VCF.
sampleIds = {} #The dictionary that relate_sample_name_to_file() will use
relate_sample_name_to_file(childFile, "child")
relate_sample_name_to_file(paternalFile, "paternal")
relate_sample_name_to_file(maternalFile, "maternal")
print(sampleIds)

#Create a dictionary that has all the variant positions of the child, for each chromosome and output a file that has variant-only positions for the child
positionDict = {} #The dictionary that filterChild() will use
filterChild(childFile)
print("Variant-only positions of the child have been written to a temporary file.")

#Output a file for each parent that has positions that occur as variant-only positions in the child
with concurrent.futures.ProcessPoolExecutor(max_workers=numberTasks) as executor:
    executor.map(filterParents, [paternalFile, maternalFile])

# Use GATK to combine all trios into one vcf and then genotype the combined trio vcf
files = ["/tmp/child_parsed.vcf.gz", "/tmp/paternal_parsed.vcf.gz", "/tmp/maternal_parsed.vcf.gz"]
tempCombinedName = "/tmp/combined.vcf.gz"
tempGenotypedName = "/tmp/genotyped.vcf.gz"
try:
    fileString = ""
    for file in files:
        fileString += f"-V {file} "
        os.system(f"gatk IndexFeatureFile -F {file}")
    # Extract fasta reference file
    os.system("unzip /fasta_references.zip -d /fasta_references")
    os.system("gzip -d /fasta_references/*.gz")
    os.system(f"gatk CombineGVCFs -R /fasta_references/Homo_sapiens_assembly38.fasta {fileString} -O {tempCombinedName}")
    print("Trio has been combined and written to a temporary file.")
    os.system(f"gatk IndexFeatureFile -F {tempCombinedName}")
    os.system(f"gatk --java-options '-Xmx2g' GenotypeGVCFs -R /fasta_references/Homo_sapiens_assembly38.fasta -V {tempCombinedName} -O {tempGenotypedName}")
    print("Trio has been joint-genotyped.")

except:
    print("Trio not combined, there was an error detected by GATK")

#Separate combined trio file by chromosome and create child scaffold
with gzip.open(tempGenotypedName, "rt") as vcf:
    outputName = f"/tmp/genotyped"
    chromosomeSet = set()
    headerChromosome = ""
    headerScaffold = ""
    fileListToBGzip = []
    for line in vcf:
        if line.startswith("##"):
            headerChromosome = headerChromosome + line
            headerScaffold = headerScaffold + line
        elif line.startswith("#CHROM"):
            headerChromosome = headerChromosome + line
            lineList = line.rstrip("\n").split("\t")
            chromIndex = lineList.index("#CHROM")
            childIndex = lineList.index(sampleIds["child"])
            paternalIndex = lineList.index(sampleIds["paternal"])
            maternalIndex = lineList.index(sampleIds["maternal"])
            newLineList = []
            for i in range(0, len(lineList)):
                if i != maternalIndex or i != paternalIndex:
                    newLineList.append(lineList[i])
            lineList = newLineList
            line = "\t".join(lineList) + "\n"
            headerScaffold = headerScaffold + line
        elif not line.startswith("#") and line.split("\t")[0] not in chromosomeSet:
            lineList = line.rstrip("\n").split("\t")
            chrom = lineList[chromIndex]
            updatedChrom = chrom[3:]
            line = line.replace(chrom, updatedChrom)
            #output metadata to output chromosome and output scaffold
            with gzip.open(f"{outputName}_{chrom}.vcf", "wb") as chromosome, gzip.open(f"{outputName}_{chrom}_scaffold.vcf", "wb") as scaffold:
                chromosome.write(headerChromosome.encode())
                chromosome.write(line.encode())
                newLineList = []
                for i in range(0, len(lineList)):
                    if i != maternalIndex or i != paternalIndex:
                        newLineList.append(lineList[i])
                lineList = newLineList
                line = "\t".join(lineList) + "\n"
                line = line.replace(chrom, updatedChrom)
                scaffold.write(headerScaffold.encode())
                scaffold.write(line.encode())
                chromosomeSet.add(chrom)
                fileListToBGzip.append(f"{outputName}_{chrom}.vcf")
                fileListToBGzip.append(f"{outputName}_{chrom}_scaffold.vcf")
        elif not line.startswith("#") and line.split("\t")[0] in chromosomeSet:
            lineList = line.rstrip("\n").split("\t")
            chrom = lineList[chromIndex]
            updatedChrom = chrom[3:]
            childGenotype = lineList[childIndex].split(":")[0]
            paternalGenotype = lineList[paternalIndex].split(":")[0]
            maternalGenotype = lineList[maternalIndex].split(":")[0]
            childAllele1 = childGenotype[0]
            childAllele2 = childGenotype[-1]
            line = line.replace(chrom, updatedChrom)
            #Output line to chromosome file
            with gzip.open(f"{outputName}_{chrom}.vcf", "ab") as chromosome:
                chromosome.write(line.encode())
            #Output line (except for parent genotypes) for trio-phased lines
            with gzip.open(f"{outputName}_{chrom}_scaffold.vcf", "ab") as scaffold:
                newLineList = []
                for i in range(0, len(lineList)):
                    if i != maternalIndex or i != paternalIndex:
                        newLineList.append(lineList[i])
                lineList = newLineList
                line = "\t".join(lineList) + "\n"
                line = line.replace(chrom, updatedChrom)
                if childAllele1 in paternalGenotype and childAllele1 not in maternalGenotype and childAllele2 in maternalGenotype:
                    phase = f"{childAllele1}|{childAllele2}"
                    line = line.replace(childGenotype, phase)
                    scaffold.write(line.encode())
                elif childAllele2 in paternalGenotype and childAllele2 not in maternalGenotype and childAllele1 in maternalGenotype:
                    phase = f"{childAllele2}|{childAllele1}"
                    line = line.replace(childGenotype, phase)
                    scaffold.write(line.encode())
                elif childAllele1 in maternalGenotype and childAllele1 not in paternalGenotype and childAllele2 in paternalGenotype:
                    phase = f"{childAllele2}|{childAllele1}"
                    line = line.replace(childGenotype, phase)
                    scaffold.write(line.encode())
                elif childAllele2 in maternalGenotype and childAllele2 not in paternalGenotype and childAllele1 in paternalGenotype:
                    phase = f"{childAllele1}|{childAllele2}"
                    line = line.replace(childGenotype, phase)
                    scaffold.write(line.encode())

for file in fileListToBGzip:
    bgzipFile(file)
    os.system(f"tabix -fp vcf {file}.gz")
    os.system(f"bcftools index {file}.gz")
    
print("Trio has been separated into chromosome files and chromosome scaffolds have been created in preparation for phasing.")

# Extract genetic maps and download haplotype references if necessary
os.system("tar -xf /shapeit4/maps/genetic_maps.b38.tar.gz -C /shapeit4/maps/")
if not os.path.exists(f"{haplotypePath}ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"):
    filesToDownload = []
    for i in range(1,23):
        filesToDownload.append(f"wget --no-check-certificate http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -P {haplotypePath}")
    with concurrent.futures.ProcessPoolExecutor(max_workers=numberTasks) as executor:
        executor.map(osSystemTask, filesToDownload)
    filesToIndex = []
    for i in range(1,23):
        filesToIndex.append(f"bcftools index {haplotypePath}ALL.chr{i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
    with concurrent.futures.ProcessPoolExecutor(max_workers=numberTasks) as executor:
        executor.map(osSystemTask, filesToIndex)

#Create a list of shapeit4 execution commands
taskList = []
for i in range(22, 0, -1):
    taskList.append(f"shapeit4 --input /tmp/genotyped_chr{i}.vcf.gz --map /shapeit4/maps/chr{i}.b38.gmap.gz --region {i} --output /tmp/phased_chr{i}_with_scaffold.vcf.gz --reference {haplotypePath}ALL.chr{i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --sequencing --scaffold /tmp/genotyped_chr{i}_scaffold.vcf.gz --seed 123456789")
#Phase with shapeit4, using concurrent.futures to phase all chromosomes at once.
with concurrent.futures.ProcessPoolExecutor(max_workers=numberTasks) as executor:
    executor.map(osSystemTask, taskList)

shapeitPositions = {}
correctlyPhased = 0
incorrectlyPhased = 0
totalVariants = 0
cannotBeDetermined = 0
allHets = 0
totalPhased = 0
header = ""
for i in range(1, 23):
    with gzip.open(f"/tmp/phased_chr{i}_with_scaffold.vcf.gz", "rt") as phasedFile:
        iteration = 0
        for line in phasedFile:
            if line.startswith("##"):
                if i == 1:
                    header = header + line
            elif line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                childIndex = lineList.index(sampleIds["child"])
                paternalIndex = lineList.index(sampleIds["paternal"])
                maternalIndex = lineList.index(sampleIds["maternal"])
                posIndex = lineList.index("POS")
                chromIndex = lineList.index("#CHROM")
                refIndex = lineList.index("REF")
                altIndex = lineList.index("ALT")
                infoIndex = lineList.index("INFO")
                if i == 1:
                    header = header + line
            else:
                totalPhased += 1
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex]
                pos = int(lineList[posIndex])
                childHaplotype = lineList[childIndex]
                paternalHaplotype = lineList[paternalIndex]
                maternalHaplotype = lineList[maternalIndex]
                childAllele1 = childHaplotype[0]
                childAllele2 = childHaplotype[-1]
                totalVariants += 1
                lineList[infoIndex] = "."
                lineList[chromIndex] = "chr" + chrom
                if chrom not in shapeitPositions:
                    shapeitPositions[chrom] = {}
                if childAllele1 in paternalHaplotype and childAllele1 not in maternalHaplotype and childAllele2 in maternalHaplotype:
                    phase = f"{childAllele1}|{childAllele2}"
                    if phase == childHaplotype:
                        correctlyPhased += 1
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    else:
                        incorrectlyPhased += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                elif childAllele2 in paternalHaplotype and childAllele2 not in maternalHaplotype and childAllele1 in maternalHaplotype:
                    phase = f"{childAllele2}|{childAllele1}"
                    if phase == childHaplotype:
                        correctlyPhased += 1
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    else:
                        incorrectlyPhased += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                elif childAllele1 in maternalHaplotype and childAllele1 not in paternalHaplotype and childAllele2 in paternalHaplotype:
                    phase = f"{childAllele2}|{childAllele1}"
                    if phase == childHaplotype:
                        correctlyPhased += 1
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    else:
                        incorrectlyPhased += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                        
                elif childAllele2 in maternalHaplotype and childAllele2 not in paternalHaplotype and childAllele1 in paternalHaplotype:
                    phase = f"{childAllele1}|{childAllele2}"
                    if phase == childHaplotype:
                        correctlyPhased += 1
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    else:
                        incorrectlyPhased += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                else:
                    line = "\t".join(lineList) + "\n"
                    shapeitPositions[chrom][pos] = line

notInShapeit = 0
for i in range(1, 23):
    with gzip.open(f"/tmp/genotyped_chr{i}.vcf.gz", "rt") as genotypeFile:
        for line in genotypeFile:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                lineList = line.rstrip("\n").split("\t")
                childIndex = lineList.index(sampleIds["child"])
                paternalIndex = lineList.index(sampleIds["paternal"])
                maternalIndex = lineList.index(sampleIds["maternal"])
                posIndex = lineList.index("POS")
                chromIndex = lineList.index("#CHROM")
                refIndex = lineList.index("REF")
                altIndex = lineList.index("ALT")
                filterIndex = lineList.index("FILTER")
                infoIndex = lineList.index("INFO")
                formatIndex = lineList.index("FORMAT")
            else:
                lineList = line.rstrip("\n").split("\t")
                chrom = lineList[chromIndex]
                pos = int(lineList[posIndex])
                childHaplotype = lineList[childIndex].split(":")[0]
                paternalGenotype = lineList[paternalIndex].split(":")[0]
                maternalGenotype = lineList[maternalIndex].split(":")[0]
                childAllele1 = childHaplotype[0]
                childAllele2 = childHaplotype[-1]
                
                if pos not in shapeitPositions[chrom]:
                    totalVariants += 1
                    lineList[filterIndex] = "."
                    lineList[infoIndex] = "."
                    lineList[formatIndex] = "GT"
                    lineList[chromIndex] = "chr" + chrom
                    lineList[paternalIndex] = paternalGenotype
                    lineList[maternalIndex] = maternalGenotype
                    if childAllele1 in paternalGenotype and childAllele1 not in maternalGenotype and childAllele2 in maternalGenotype:
                        phase = f"{childAllele1}|{childAllele2}"
                        notInShapeit += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    elif childAllele2 in paternalGenotype and childAllele2 not in maternalGenotype and childAllele1 in maternalGenotype:
                        phase = f"{childAllele2}|{childAllele1}"
                        notInShapeit += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    elif childAllele1 in maternalGenotype and childAllele1 not in paternalGenotype and childAllele2 in paternalGenotype:
                        phase = f"{childAllele2}|{childAllele1}"
                        notInShapeit += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    elif childAllele2 in maternalGenotype and childAllele2 not in paternalGenotype and childAllele1 in paternalGenotype:
                        phase = f"{childAllele1}|{childAllele2}"
                        notInShapeit += 1
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line
                    elif childAllele1 == childAllele2 and (childAllele1 in paternalGenotype and childAllele1 in maternalGenotype):
                        phase = f"{childAllele1}|{childAllele2}"
                        lineList[childIndex] = phase
                        line = "\t".join(lineList) + "\n"
                        shapeitPositions[chrom][pos] = line

print(sampleIds)
print(f"There were {correctlyPhased} ({(correctlyPhased / (correctlyPhased + incorrectlyPhased)) * 100:.2f}%) correctly phased haplotypes, and {incorrectlyPhased} ({(incorrectlyPhased / (correctlyPhased + incorrectlyPhased)) * 100:.2f}%) incorrectly phased haplotypes. {notInShapeit} variants were not phased by shapeit but were phasable and will be included in final output.")
print(f"There were {totalVariants} total variants phased.\n")

with gzip.open(outputFile.replace(".gz", ""), "wb") as output:
    output.write(header.encode())
    for chrom, posDict in sorted(shapeitPositions.items()):
        for pos, line in sorted(posDict.items()):
            output.write(line.encode())

os.system(f"zcat {outputFile.replace('.gz', '')} | bgzip -f > {outputFile}")
os.system(f"tabix -fp vcf {outputFile}")
os.system(f"bcftools index {outputFile}")
os.system(f"rm {outputFile.replace('.gz', '')}")
print(f"Phased output file written as {outputFile}")

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours) {char}')