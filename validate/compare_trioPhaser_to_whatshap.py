import gzip
import argparse

# Argparse Information
parser = argparse.ArgumentParser(description='Compares the trioPhaser phased \
                                output to the whatshap phased output')
parser.add_argument('trioPhaser_output_file', 
                    help='The location and name of the trioPhaser output file.')
parser.add_argument('whatshap_output_file', 
                    help='The location and name of the whatshap output file.')
parser.add_argument('output_files_path', 
                    help = 'The path where the output files should be stored.\
                        There is an output for congrugent positions, and output\
                        for positions phased by trioPhaser but not whatshap, a file \
                        for positions phased incorrectly by whatshap')

args = parser.parse_args()

# Create variables of each argument from argparse
trioPhaser_file = args.trioPhaser_output_file
whatshap_file = args.whatshap_output_file
output_path = args.output_files_path

# Functions
def get_header_indexes(header_line_list):
    chrom_index = header_line_list.index("#CHROM")
    pos_index = header_line_list.index("POS")
    ref_index = header_line_list.index("REF")
    alt_index = header_line_list.index("ALT")
    qual_index = header_line_list.index("QUAL")
    return(chrom_index, pos_index, ref_index, alt_index, qual_index)

def get_header_info(line_list, chrom_index, pos_index, ref_index, alt_index, qual_index):
    chrom = line_list[chrom_index]
    pos = line_list[pos_index]
    ref = line_list[ref_index]
    alt = line_list[alt_index]
    qual = line_list[qual_index]
    return(chrom, pos, ref, alt, qual)

def nucleotide_based_haplotype(input_list, ref_allele, alt_allele_list):
    """ Get a nucleotide-based haplotype (i.e. A|G) instead of a numeric-based 
        haplotype (i.e. 0|1).
    """
    temp_list = []
    for i in input_list:
        if i == "0":
            temp_list.append(ref_allele)
        else:
            temp_list.append(alt_allele_list[int(i) - 1])
    temp_haplotype = "|".join(temp_list)
    return(temp_haplotype)

# Create a dictionary of all positions phased by whatshap with a 
# quality >= 30.
haplotype_dict = {}
whatshap_total_variants = 0
with open(whatshap_file, "rt") as whatshap:
    for line in whatshap:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            if "21228" in line:
                sample_index = line_list.index("21228")
            elif "38536" in line:
                sample_index = line_list.index("38536")
        else:
            line_list = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, qual = get_header_info(line_list, chrom_index, 
                                                        pos_index, ref_index, 
                                                        alt_index, qual_index)
            sample_column = line_list[sample_index]
            haplotype = sample_column[0:3]
            sample_columnList = sample_column.split(":")
            if chrom[3:].isnumeric() and int(chrom[3:]) in range(1, 23):
                if chrom not in haplotype_dict \
                    and "|" in haplotype \
                    and "." not in haplotype \
                    and qual != "." \
                    and float(qual) >= 30:
                    whatshap_total_variants += 1
                    haplotype_dict[chrom] = {pos: [ref, alt, haplotype]}
                elif chrom in haplotype_dict \
                    and "|" in haplotype \
                    and "." not in haplotype \
                    and qual != "." \
                    and float(qual) >= 30:
                    whatshap_total_variants += 1
                    haplotype_dict[chrom][pos] = [ref, alt, haplotype]

# Compare the positions phased by trioPhaser with the positions phased by whatshap's 
# Long Ranger.

phased_correctly = 0
phased_incorrectly = 0
phased_but_not_in_whatshap = 0
different_calls = 0
mendelian_phased = 0
mendelian_phased_wrong = 0
mendelian_phased_total = 0
trio_phaser_total_phased = 0
mendelian_phasable = {}
with gzip.open(trioPhaser_file, "rt") as trio_phaser:
    for line in trio_phaser:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            if "21228" in line:
                sample_index = line_list.index("21228")
                paternal_index = line_list.index("21230")
                maternal_index = line_list.index("21229")
            elif "38536" in line:
                sample_index = line_list.index("38536")
                paternal_index = line_list.index("38535")
                maternal_index = line_list.index("38534")
        else:
            line_list = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, qual = get_header_info(line_list, chrom_index, 
                                                        pos_index, ref_index, 
                                                        alt_index, qual_index)
            alt_list = alt.split(",")
            haplotype = line_list[sample_index]
            haplotype_list = haplotype.split("|")
            paternal_haplotype = line_list[paternal_index]
            if "|" in paternal_haplotype:
                paternal_haplotype_list = paternal_haplotype.split("|")
            else:
                paternal_haplotype_list = paternal_haplotype.split("/")
            maternal_haplotype = line_list[maternal_index]
            if "|" in maternal_haplotype:
                maternal_haplotype_list = maternal_haplotype.split("|")
            else:
                maternal_haplotype_list = maternal_haplotype.split("/")

            # This accounts for multiallelic positions and creates a bi-allelic
            new_haplotype_child = nucleotide_based_haplotype(haplotype_list, ref, alt_list)
            new_haplotype_paternal = nucleotide_based_haplotype(paternal_haplotype_list, ref, alt_list)
            new_haplotype_maternal = nucleotide_based_haplotype(maternal_haplotype_list, ref, alt_list)

            child_allele_1 = new_haplotype_child.split("|")[0]
            child_allele_2 = new_haplotype_child.split("|")[-1]

            if child_allele_1 in new_haplotype_paternal \
                and child_allele_1 not in new_haplotype_maternal \
                and child_allele_2 in new_haplotype_maternal:
                phase = f"{child_allele_1}|{child_allele_2}"
            elif child_allele_2 in new_haplotype_paternal \
                and child_allele_2 not in new_haplotype_maternal \
                and child_allele_1 in new_haplotype_maternal:
                phase = f"{child_allele_2}|{child_allele_1}"
            elif child_allele_1 in new_haplotype_maternal \
                and child_allele_1 not in new_haplotype_paternal \
                and child_allele_2 in new_haplotype_paternal:
                phase = f"{child_allele_2}|{child_allele_1}"
            elif child_allele_2 in new_haplotype_maternal \
                and child_allele_2 not in new_haplotype_paternal \
                and child_allele_1 in new_haplotype_paternal:
                phase = f"{child_allele_1}|{child_allele_2}"
            elif child_allele_1 == child_allele_2 \
                and (child_allele_1 in new_haplotype_maternal \
                and child_allele_1 in new_haplotype_paternal):
                phase = f"{child_allele_1}|{child_allele_2}"
            else:
                phase = "."

            if phase != "." and child_allele_1 != child_allele_2:
                mendelian_phased_total += 1
                trio_phaser_total_phased += 1
                if chrom not in mendelian_phasable:
                    mendelian_phasable[chrom] = set(pos)
                else:
                    mendelian_phasable[chrom].add(pos)
            elif phase == "." and child_allele_1 != child_allele_2:
                trio_phaser_total_phased += 1

phased_incorrectly_by_whatshap = {}
with gzip.open(trioPhaser_file, "rt") as trio_phaser, \
    open(f"{output_path}issues.tsv", 'wt') as output, \
    gzip.open(f"{output_path}phased_by_trio_phaser_but_not_whatshap.vcf.gz", "wb") as unique_to_trio_phaser, \
    gzip.open(f"{output_path}phased_congruently.vcf.gz", "wb") as phased_out:
    for line in trio_phaser:
        if line.startswith("##"):
            unique_to_trio_phaser.write(line.encode())
            phased_out.write(line.encode())
        elif line.startswith("#CHROM"):

            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            if "21228" in line:
                sample_index = line_list.index("21228")
            elif "38536" in line:
                sample_index = line_list.index("38536")
            output.write("chrom\tpos\twhatshap\ttrio_phaser\n")
            unique_to_trio_phaser.write(line.encode())
            phased_out.write(line.encode())
        else:
            line_list = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, qual = get_header_info(line_list, chrom_index, 
                                                        pos_index, ref_index, 
                                                        alt_index, qual_index)
            alt_list = alt.split(",")
            haplotype = line_list[sample_index]
            haplotype_list = haplotype.split("|")
            if pos in haplotype_dict[chrom]:
                ref_whatshap = haplotype_dict[chrom][pos][0]
                alt_whatshap = haplotype_dict[chrom][pos][1]
                alt_whatshap_list = alt_whatshap.split(",")
                haplotype_whatshap = haplotype_dict[chrom][pos][-1]
                haplotype_whatshap_list = haplotype_whatshap.split("|")
                
                # Get a nucleotide-based haplotype instead of a numeric-based 
                # haplotype for the whatshap data.
                new_haplotype_whatshap = nucleotide_based_haplotype(haplotype_whatshap_list, ref_whatshap, alt_whatshap_list)

                # Get a nucleotide-based haplotype instead of a numeric-based 
                # haplotype for the trio_phaser data.
                new_haplotype = nucleotide_based_haplotype(haplotype_list, ref, alt_list)
            
                if new_haplotype_whatshap == new_haplotype \
                    and pos in mendelian_phasable[chrom]:
                    phased_correctly += 1
                    mendelian_phased += 1
                    phased_out.write(line.encode())
                elif new_haplotype_whatshap == new_haplotype \
                    and pos not in mendelian_phasable[chrom]:
                    phased_correctly += 1
                    phased_out.write(line.encode())
                elif f"{new_haplotype_whatshap.split('|')[-1]}|{new_haplotype_whatshap.split('|')[0]}" == new_haplotype \
                    and new_haplotype != new_haplotype_whatshap and pos in mendelian_phasable[chrom]:
                    phased_incorrectly += 1
                    mendelian_phased_wrong += 1
                    if chrom not in phased_incorrectly_by_whatshap:
                        phased_incorrectly_by_whatshap[chrom] = [pos]
                    else:
                        phased_incorrectly_by_whatshap[chrom].append(pos)
                    output.write(f"{chrom}\t{pos}\t{new_haplotype_whatshap}\t{new_haplotype}\n")
                elif f"{new_haplotype_whatshap.split('|')[-1]}|{new_haplotype_whatshap.split('|')[0]}" == new_haplotype \
                    and new_haplotype != new_haplotype_whatshap and pos not in mendelian_phasable[chrom]:
                    phased_incorrectly += 1
                    #output.write(f"{chrom}\t{pos}\t{haplotype_dict[chrom][pos]}\t{[ref, alt, haplotype]}\n")
                elif f"{new_haplotype_whatshap.split('|')[-1]}|{new_haplotype_whatshap.split('|')[0]}" != new_haplotype \
                    and new_haplotype_whatshap != new_haplotype:
                    different_calls += 1
            elif pos not in haplotype_dict[chrom]:
                if haplotype_list[0] != haplotype_list[-1]:
                    unique_to_trio_phaser.write(line.encode())
                    phased_but_not_in_whatshap += 1

with open(whatshap_file, "rt") as input, \
    gzip.open(f"{output_path}mendelian_phased_incorrectly_by_whatshap.vcf.gz", "wb") as mendelian_out:
     for line in input:
        if line.startswith("##"):
            mendelian_out.write(line.encode())
        elif line.startswith("#CHROM"):
            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            mendelian_out.write(line.encode())
        else:
            line_list = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, qual = get_header_info(line_list, chrom_index, 
                                                        pos_index, ref_index, 
                                                        alt_index, qual_index)
            if chrom in phased_incorrectly_by_whatshap \
                and pos in phased_incorrectly_by_whatshap[chrom]:
                mendelian_out.write(line.encode())


print(f"There were {whatshap_total_variants} phased variants by whatshap.")
print(f"There were {trio_phaser_total_phased} phased variants by trio_phaser.")
print(f"There were {phased_correctly + phased_incorrectly} phased variants that were able to be compared between trio_phaser and whatshap (i.e. the ref and alt(s) were called the same at a given position).")
print(f"There were {phased_correctly} ({(phased_correctly / (phased_correctly + phased_incorrectly)) * 100}%) variants phased the same between whatshap and trio_phaser")
print(f"There were {phased_but_not_in_whatshap} variants phased through trio_phaser that were not phased in whatshap")
print(f'There were {different_calls} variants that were genotyped differently between whatshap and GATK')
print(f'{mendelian_phased} {100 - ((mendelian_phased_wrong / mendelian_phased_total) * 100)}% of the Mendelian variants were phased correctly by whatshap.')
print(f'There were {mendelian_phased_wrong} Mendelian variants that were incorrectly phased by whatshap.')