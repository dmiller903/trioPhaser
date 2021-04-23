import gzip

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

# Create a dictionary of all positions phased by 10X Long Ranger with a 
# quality >= 20.
haplotype_dict = {}
phase_set = {}
tenX_total_variants = 0
with gzip.open("son_GRCh38_longRanger.vcf.gz", "rt") as long_ranger:
    for line in long_ranger:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            sample_index = line_list.index("21228")
        else:
            line_list = line.rstrip("\n").split("\t")
            chrom = line_list[chrom_index]
            pos = line_list[pos_index]
            ref = line_list[ref_index]
            alt = line_list[alt_index]
            qual = line_list[qual_index]
            sample_column = line_list[sample_index]
            haplotype = sample_column[0:3]
            sample_columnList = sample_column.split(":")
            phase_set_value = sample_columnList[-3]
            if chrom[3:].isnumeric() and int(chrom[3:]) in range(1, 23):
                if chrom not in haplotype_dict \
                    and "|" in haplotype \
                    and qual != "." \
                    and float(qual) >= 20:
                    tenX_total_variants += 1
                    haplotype_dict[chrom] = {pos: [ref, alt, haplotype]}
                    phase_set[chrom] = {pos: phase_set_value}
                elif chrom in haplotype_dict \
                    and "|" in haplotype \
                    and qual != "." \
                    and float(qual) >= 20:
                    tenX_total_variants += 1
                    haplotype_dict[chrom][pos] = [ref, alt, haplotype]
                    phase_set[chrom][pos] = phase_set_value

# Compare the positions phased by trioPhaser with the positions phased by 10X's 
# Long Ranger.

phased_correctly = 0
phased_incorrectly = 0
phased_but_not_in_10X = 0
different_calls = 0
mendelian_phased = 0
mendelian_phased_wrong = 0
mendelian_phased_total = 0
phase_set_number_wrong = {}
number_in_phase_set = {}
trio_phaser_total_phased = 0
mendelian_phasable = {}
with gzip.open("giab_phased.vcf.gz", "rt") as trio_phaser:
    for line in trio_phaser:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            sample_index = line_list.index("21228")
            paternal_index = line_list.index("21230")
            maternal_index = line_list.index("21229")
        else:
            trio_phaser_total_phased += 1
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

            if phase != ".":
                mendelian_phased_total += 1
                if chrom not in mendelian_phasable:
                    mendelian_phasable[chrom] = set(pos)
                else:
                    mendelian_phasable[chrom].add(pos)
            
            if pos in haplotype_dict[chrom]:
                ref_10X = haplotype_dict[chrom][pos][0]
                alt_10X = haplotype_dict[chrom][pos][1]
                alt_10X_list = alt_10X.split(",")
                haplotype_10X = haplotype_dict[chrom][pos][-1]
                haplotype_10X_list = haplotype_10X.split("|")
                
                # Get a nucleotide-based haplotype instead of a numeric-based 
                # haplotype for the 10X data.
                new_haplotype_10X = nucleotide_based_haplotype(haplotype_10X_list, ref_10X, alt_10X_list)

                # Get a nucleotide-based haplotype instead of a numeric-based 
                # haplotype for the trio_phaser data.
                new_haplotype = nucleotide_based_haplotype(haplotype_list, ref, alt_list)
                
                phase_set_value = phase_set[chrom][pos]
                if phase != ".":
                    if new_haplotype == phase \
                        and new_haplotype != new_haplotype_10X \
                        and f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype:
                        if phase_set_value not in phase_set_number_wrong:
                            phase_set_number_wrong[phase_set_value] = 1
                            number_in_phase_set[phase_set_value] = 1
                        elif phase_set_value in phase_set_number_wrong:
                            phase_set_number_wrong[phase_set_value] += 1
                            number_in_phase_set[phase_set_value] += 1
                    elif new_haplotype == new_haplotype_10X \
                        and new_haplotype == phase:
                        if phase_set_value not in number_in_phase_set:
                            number_in_phase_set[phase_set_value] = 1
                        elif phase_set_value in number_in_phase_set:
                            number_in_phase_set[phase_set_value] += 1

phased_incorrectly_by_10X = {}
with gzip.open("giab_phased.vcf.gz", "rt") as trio_phaser, \
    open("issues.tsv", 'wt') as output, \
    gzip.open("phased_by_trio_phaser_but_not_10X.vcf.gz", "wb") as unique_to_trio_phaser, \
    gzip.open("phased_congruently.vcf.gz", "wb") as phased_out:
    for line in trio_phaser:
        if line.startswith("##"):
            unique_to_trio_phaser.write(line.encode())
            phased_out.write(line.encode())
        elif line.startswith("#CHROM"):

            line_list = line.rstrip("\n").split("\t")
            chrom_index, pos_index, ref_index, \
            alt_index, qual_index = get_header_indexes(line_list)
            sample_index = line_list.index("21228")
            output.write("chrom\tpos\t10x\tphase_set_value\ttrio_phaser\n")
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
                ref_10X = haplotype_dict[chrom][pos][0]
                alt_10X = haplotype_dict[chrom][pos][1]
                alt_10X_list = alt_10X.split(",")
                haplotype_10X = haplotype_dict[chrom][pos][-1]
                haplotype_10X_list = haplotype_10X.split("|")
                
                # Get a nucleotide-based haplotype instead of a numeric-based 
                # haplotype for the 10X data.
                new_haplotype_10X = nucleotide_based_haplotype(haplotype_10X_list, ref_10X, alt_10X_list)

                # Get a nucleotide-based haplotype instead of a numeric-based 
                # haplotype for the trio_phaser data.
                new_haplotype = nucleotide_based_haplotype(haplotype_list, ref, alt_list)
            
                phase_set_value = phase_set[chrom][pos]
                if phase_set_value in phase_set_number_wrong:
                    if (phase_set_number_wrong[phase_set_value] / number_in_phase_set[phase_set_value]) > 0.51:
                        if f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype \
                            and pos in mendelian_phasable[chrom]:
                            phased_correctly += 1
                            mendelian_phased += 1
                            phased_out.write(line.encode())
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype \
                            and pos not in mendelian_phasable[chrom]:
                            phased_correctly += 1
                            phased_out.write(line.encode())
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" != new_haplotype \
                            and new_haplotype == new_haplotype_10X \
                            and pos in mendelian_phasable[chrom]:
                            phased_incorrectly += 1
                            mendelian_phased_wrong += 1
                            if chrom not in phased_incorrectly_by_10X:
                                phased_incorrectly_by_10X[chrom] = [pos]
                            else:
                                phased_incorrectly_by_10X[chrom].append(pos)
                            output.write(f"{chrom}\t{pos}\t{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}\t{phase_set[chrom][pos]}\t{new_haplotype}\n")
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" != new_haplotype \
                            and new_haplotype == new_haplotype_10X \
                            and pos not in mendelian_phasable[chrom]:
                            phased_incorrectly += 1
                            #output.write(f"{chrom}\t{pos}\t{haplotype_dict[chrom][pos]}\t{phase_set[chrom][pos]}\t{[ref, alt, haplotype]}\n")
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" != new_haplotype \
                            and new_haplotype_10X != new_haplotype:
                            different_calls += 1
                    else:
                        if new_haplotype_10X == new_haplotype \
                            and pos in mendelian_phasable[chrom]:
                            phased_correctly += 1
                            mendelian_phased += 1
                            phased_out.write(line.encode())
                        elif new_haplotype_10X == new_haplotype \
                            and pos not in mendelian_phasable[chrom]:
                            phased_correctly += 1
                            phased_out.write(line.encode())
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype \
                            and new_haplotype != new_haplotype_10X and pos in mendelian_phasable[chrom]:
                            phased_incorrectly += 1
                            mendelian_phased_wrong += 1
                            if chrom not in phased_incorrectly_by_10X:
                                phased_incorrectly_by_10X[chrom] = [pos]
                            else:
                                phased_incorrectly_by_10X[chrom].append(pos)
                            output.write(f"{chrom}\t{pos}\t{new_haplotype_10X}\t{phase_set[chrom][pos]}\t{new_haplotype}\n")
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype \
                            and new_haplotype != new_haplotype_10X \
                            and pos not in mendelian_phasable[chrom]:
                            phased_incorrectly += 1
                            #output.write(f"{chrom}\t{pos}\t{haplotype_dict[chrom][pos]}\t{phase_set[chrom][pos]}\t{[ref, alt, haplotype]}\n")
                        elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" != new_haplotype \
                            and new_haplotype_10X != new_haplotype:
                            different_calls += 1
                else:
                    if new_haplotype_10X == new_haplotype \
                        and pos in mendelian_phasable[chrom]:
                        phased_correctly += 1
                        mendelian_phased += 1
                        phased_out.write(line.encode())
                    elif new_haplotype_10X == new_haplotype \
                        and pos not in mendelian_phasable[chrom]:
                        phased_correctly += 1
                        phased_out.write(line.encode())
                    elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype \
                        and new_haplotype != new_haplotype_10X and pos in mendelian_phasable[chrom]:
                        phased_incorrectly += 1
                        mendelian_phased_wrong += 1
                        if chrom not in phased_incorrectly_by_10X:
                            phased_incorrectly_by_10X[chrom] = [pos]
                        else:
                            phased_incorrectly_by_10X[chrom].append(pos)
                        output.write(f"{chrom}\t{pos}\t{new_haplotype_10X}\t{phase_set[chrom][pos]}\t{new_haplotype}\n")
                    elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" == new_haplotype \
                        and new_haplotype != new_haplotype_10X and pos not in mendelian_phasable[chrom]:
                        phased_incorrectly += 1
                        #output.write(f"{chrom}\t{pos}\t{haplotype_dict[chrom][pos]}\t{phase_set[chrom][pos]}\t{[ref, alt, haplotype]}\n")
                    elif f"{new_haplotype_10X.split('|')[-1]}|{new_haplotype_10X.split('|')[0]}" != new_haplotype \
                        and new_haplotype_10X != new_haplotype:
                        different_calls += 1
            else:
                unique_to_trio_phaser.write(line.encode())
                phased_but_not_in_10X += 1

with gzip.open("son_GRCh38_longRanger.vcf.gz", "rt") as input, \
    gzip.open("mendelian_phased_incorrectly_by_10X.vcf.gz", "wb") as mendelian_out:
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
            if chrom in phased_incorrectly_by_10X \
                and pos in phased_incorrectly_by_10X[chrom]:
                mendelian_out.write(line.encode())


print(f"There were {tenX_total_variants} phased variants by 10X.")
print(f"There were {trio_phaser_total_phased} phased variants by trio_phaser.")
print(f"There were {phased_correctly + phased_incorrectly} phased variants that were able to be compared between trio_phaser and 10X (i.e. the ref and alt(s) were called the same at a given position).")
print(f"There were {phased_correctly} ({(phased_correctly / (phased_correctly + phased_incorrectly)) * 100}%) variants phased the same between 10X and trio_phaser")
print(f"There were {phased_but_not_in_10X} variants phased through trio_phaser that were not phased in 10X")
print(f'There were {different_calls} variants that were genotyped differently between 10X and GATK')
print(f'{mendelian_phased} {100 - ((mendelian_phased_wrong / mendelian_phased_total) * 100)}% of the Mendelian variants were phased correctly by 10X.')
print(f'There were {mendelian_phased_wrong} Mendelian variants that were incorrectly phased by 10X.')