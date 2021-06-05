import re
import glob
import argparse

# Argparse Information
parser = argparse.ArgumentParser(description='Averages the summary statistics\
    across all Neuroblastoma samples')
parser.add_argument('path_to_family_directories', 
                    help='The path were all trio/family directories are stored')
parser.add_argument('name_of_output', 
                    help='The name of the output file. The output file\
                        contains all the information output by the \
                        "trio_phaser.py" script as a single tsv file.')

args = parser.parse_args()

# Create variables of each argument from argparse
input_path = args.path_to_family_directories
output_path = args.name_of_output

stat_summary = {"initial variants": 0, "correctly phased by SHAPEIT4": 0, 
                "incorrectly phased by SHAPEIT4": 0, 
                "phaseable but not phased by SHAPEIT4": 0, 
                "phased by SHAPEIT4 but could not be verified": 0, 
                "total phased": 0, "phased by SHAPEIT4": 0, 
                "time elapsed (hours)": 0, 
                "het phaseable but not phased by SHAPEIT4": 0}
with open(output_path, "wt") as output:
    output.write("family_id\tvalue\tstage\n")
    for file in glob.glob(f"{input_path}FM_*/trio_phaser_*.out"):
        with open(file) as inputFile:
            for line in inputFile:
                # Get family ID
                if "Positions in FM_" in line:
                    family_id = re.findall(r"Positions in (FM_\w+)\/", line)[0]
                # This outputs the number of initial variants for current 
                # sample and adds those variants to the overall tally in the 
                # "stat_summary" dict.
                if "Traversal complete. Processed" in line:
                    positions = re.findall(r"Traversal complete\. Processed (\d+) total", line)[0]
                    initial_variants = int(positions)
                if "Successfully wrote index to /tmp/child_parsed.vcf.gz.tbi" in line:
                    output.write(f"{family_id}\t{initial_variants}\tinitial variants\n")
                    stat_summary["initial variants"] += initial_variants

                # This outputs the number of correctly and incorrectly phased
                # variants by SHAPEIT4. It adds the number of correctly and 
                # incorrectly phased variants to the overall tally in the 
                # "stat_summary" dict.
                if "correctly phased haplotypes, and " in line:
                    group_list = re.findall(r"There were (\d+) \([\d\.]+%\) correctly phased haplotypes, and (\d+) \([\d|\.]+%\) incorrectly phased haplotypes\.", line)[0]
                    correctly_phased = int(group_list[0])
                    incorrectly_phased = int(group_list[1])
                    output.write(f"{family_id}\t{correctly_phased}\tcorrectly phased by SHAPEIT4\n")
                    stat_summary["correctly phased by SHAPEIT4"] += correctly_phased
                    output.write(f"{family_id}\t{incorrectly_phased}\tincorrectly phased by SHAPEIT4\n")
                    stat_summary["incorrectly phased by SHAPEIT4"] += incorrectly_phased
                
                # This outputs the number of phaseable (by Mendelian inhertiance)
                # variants that were not phased by SHAPEIT4 and adds those 
                # variants to the overall tally in the "stat_summary" dict.
                if "variants were not phased by shapeit but were phaseable" in line:
                    phaseable = int(re.findall(r"(\d+) variants were not phased by shapeit but were phaseable and will be included in final output\.", line)[0])
                    output.write(f"{family_id}\t{phaseable}\tphased by Mendelian inheritance\n")
                    stat_summary["phaseable but not phased by SHAPEIT4"] += phaseable

                # This outputs the number of phaseable het (by Mendelian inhertiance)
                # variants that were not phased by SHAPEIT4 and adds those 
                # variants to the overall tally in the "stat_summary" dict.
                if "variants not phased by shapeit but were phaseable," in line:
                    het = int(re.findall(r"variants not phased by shapeit but were phaseable, (\d+)", line)[0])
                    output.write(f"{family_id}\t{het}\thet phased by Mendelian inheritance\n")
                    stat_summary["het phaseable but not phased by SHAPEIT4"] += het

                # This outputs the number of variants that were phased by 
                # SHAPEIT4 but were not able to be verified with Mendelian 
                # inheritance and adds those variants to the overall tally in 
                # the "stat_summary" dict.
                if "could not be verified by using Mendelian inheritance alone." in line:
                    not_verified = int(re.findall(r"(\d+) phased variants \(as phased by SHAPEIT4\) could not be verified by using Mendelian inheritance alone\.", line)[0])
                    output.write(f"{family_id}\t{not_verified}\tphased by SHAPEIT4 but could not be verified\n")
                    stat_summary["phased by SHAPEIT4 but could not be verified"] += not_verified

                # This outputs the total number of phased variants by trioPhaser
                # and the total number of variants that were phased with 
                # SHAPEIT4. It then adds the total number of phased variants by 
                # trioPhaser and phased variants by SHAPEIT4 to the overall 
                # tally in the "stat_summary" dict.
                if "total variants phased." in line:
                    total_phased = int(re.findall(r"There were (\d+) total variants phased\.", line)[0])
                    phased_by_shapeit4 = total_phased - phaseable
                    output.write(f"{family_id}\t{total_phased}\ttotal phased with trioPhaser\n")
                    output.write(f"{family_id}\t{phased_by_shapeit4}\tphased by SHAPEIT4\n")
                    stat_summary["total phased"] += total_phased
                    stat_summary["phased by SHAPEIT4"] += phased_by_shapeit4
                
                # This outputs the number of hours it took for trioPhaser to
                # phase a sample. It then adds the number of hours to the
                # overall tally in the "stat_summary" dict.
                if "Done. Time elapsed:" in line:
                    hours = float(re.findall(r"Done\. Time elapsed: [\d\.]+ minutes \(([\d\.]+) hours\)", line)[0])
                    output.write(f"{family_id}\t{hours}\ttime elapsed (hours)\n")
                    stat_summary["time elapsed (hours)"] += hours
                
# Print an average for the value of each key in stat_summary
for key, value in stat_summary.items():
    print(f"{key}: {value / 50}")