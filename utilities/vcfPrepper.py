import os
import subprocess
import sys

# Get the directory of the current script
current_dir = os.path.dirname(os.path.abspath(__file__))

# Paths for VCF files and output
vcfFileFolderPath = current_dir  # Path where VCF folder is present
outputFolderPath = os.path.join(current_dir, "output")  # Output folder in the same directory as the script
snpEffPath = os.path.join(vcfFileFolderPath, "utilities", "SnpEff")  # Path of SnpEff folder
databasePath = os.path.join(vcfFileFolderPath, "gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz")  # Path of the database file

genome = 'GRCh38.p13'

snpEffJarPath = os.path.join(snpEffPath, 'snpEff.jar')
snpSiftJarPath = os.path.join(snpEffPath, 'snpSift.jar')

# Change command based on the operating system
if sys.platform == "Windows":
    MOVE = "move"
    MKDIR = "mkdir"
else:
    MOVE = "mv"
    MKDIR = "md"

os.chdir(snpEffPath)

# Ensure output folder exists
os.makedirs(outputFolderPath, exist_ok=True)

# Process each VCF file in the specified folder
for vcfFile in os.listdir(vcfFileFolderPath):
    if vcfFile.endswith(".vcf") or vcfFile.endswith("g.vcf"):  # Ensure only VCF files are processed
        print('Processing file: ' + vcfFile)
        filePath = os.path.join(vcfFileFolderPath, vcfFile)  # Use os.path.join for compatibility
        outputPath = os.path.join(outputFolderPath, vcfFile)  # Output path for processed file

        # Prepare command for SnpEff
        snpEffCmd = ['java', '-jar', snpEffJarPath, '-o', 'vcf', genome, filePath]

        # Run SnpEff
        intVCFPath = os.path.join(snpEffPath, 'intVCF.vcf')
        with open(intVCFPath, 'w') as intVCF:
            subprocess.run(snpEffCmd, stdout=intVCF)

        # Process the intermediate VCF file
        changedAFPath = os.path.join(snpEffPath, 'changedAF.vcf')
        with open(intVCFPath, "r") as orig, open(changedAFPath, 'w') as updated:
            for line in orig:
                newString = line.replace(';AF=', ';CHANGEDAF=')
                updated.write(newString)

        # Annotate the VCF using SnpSift
        snpSiftCmd = ['java', '-jar', snpSiftJarPath, 'annotate', databasePath, '-info', 'AF', changedAFPath]

        annotatedVCFPath = os.path.join(snpEffPath, 'annotatedVCF.vcf')

        with open(annotatedVCFPath, 'w') as annotatedFile:
            subprocess.run(snpSiftCmd, stdout=annotatedFile)

        # Final processing and renaming of files
        with open(annotatedVCFPath, "rt") as orig, open(outputPath, "wt") as updated:
            for line in orig:
                newLine = line.replace(';AF=', ';POPAF=')
                newLine = newLine.replace(';CHANGEDAF=', ';AF=')
                updated.write(newLine)

        # Clean up temporary files
        os.remove(intVCFPath)
        os.remove(annotatedVCFPath)
        os.remove(changedAFPath)

print("Processing completed.")
