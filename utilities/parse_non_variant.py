import json
import vcf
from tqdm import tqdm
from ref import extract_chromosome_sequence


def normalize_chromosome(chrom):
    """
    Normalizes chromosome names to the `Chr` format.

    Args:
        chrom (str): Original chromosome name.

    Returns:
        str: Normalized chromosome name (e.g., `Chr1`, `ChrM`, `ChrX`).
    """
    if chrom.startswith("chr"):
        chrom_tag = chrom[3:]  # Remove "chr" prefix
    else:
        chrom_tag = chrom

    # Handle special cases for mitochondrion
    if chrom_tag.lower() in ["m", "mitochondrion"]:
        return "ChrM"
    elif chrom_tag.upper() == "X":
        return "ChrX"
    elif chrom_tag.upper() == "Y":
        return "ChrY"
    else:
        return f"Chr{chrom_tag}"


def determine_allelic_state(chrom, gt):
    """
    Determines the allelic state based on the chromosome type and genotype.

    Args:
        chrom (str): Chromosome name.
        gt (str): Genotype (e.g., "0/0" or "0").

    Returns:
        str: Allelic state (`'homoplasmic'`, `'homozygous'`, `'hemizygous'`).
    """
    if chrom == "ChrM":  # For mitochondrial chromosomes
        if gt in ["0/0", "0"]:
            return "homoplasmic"
    elif gt == "0/0":  # For autosomal chromosomes
        return "homozygous"
    elif gt == "0":  # For autosomal chromosomes
        return "hemizygous"

    return "unknown"  # In case the genotype doesn't match expected patterns


def parse_gvcf(gvcf_filename, build_version):
    """
    Parses the GVCF file and identifies contiguous studied regions.

    Args:
        gvcf_filename (str): Path to the GVCF file.
        build_version (int): Build version of the reference genome (37 or 38).

    Returns:
        tuple: A list of regions (BED format) and a list of detailed records.
    """
    studied_regions = []  # List to store contiguous regions as dictionaries
    detailed_blocks = []  # To store detailed records

    # Read the GVCF file
    vcf_reader = vcf.Reader(filename=gvcf_filename)

    current_chromosome = None
    current_start = None
    current_end = None

    for record in vcf_reader:
        chrom = normalize_chromosome(record.CHROM)
        pos = record.POS  # 1-based position
        end = record.INFO.get('END', pos)  # Use 'END' field if available, else default to 'POS'

        # Initialize or extend the current region
        if current_start is None:
            current_chromosome = chrom
            current_start = pos
            current_end = end
        elif chrom == current_chromosome and pos <= current_end + 1:
            # Extend the region
            current_end = max(current_end, end)
        else:
            # Finalize the current region and start a new one
            studied_regions.append({
                "CHROM": current_chromosome,
                "START": current_start - 1,  # Convert to 0-based for BED
                "END": current_end  # Keep as 1-based
            })
            current_chromosome = chrom
            current_start = pos
            current_end = end

        # Add detailed block metadata, including allelicState
        gt = "0" if chrom == "ChrM" else record.samples[0].data.GT
        allelic_state = determine_allelic_state(chrom, gt)

        detailed_blocks.append({
            "CHROM": chrom,
            "POS": pos,
            "END": end,
            "REF_ALLELE": record.REF,
            "FILTER": record.FILTER,
            "GT": gt,
            "allelicState": allelic_state  # Add allelicState to the block
        })

    # Finalize the last region
    if current_start is not None:
        studied_regions.append({
            "CHROM": current_chromosome,
            "START": current_start - 1,  # Convert to 0-based for BED
            "END": current_end  # Keep as 1-based
        })

    return studied_regions, detailed_blocks


def generate_bed_json(studied_regions, build_version):
    """
    Generates BED JSON for contiguous studied regions.

    Args:
        studied_regions (dict): Studied regions per chromosome.
        build_version (int): Build version of the reference genome (37 or 38).

    Returns:
        list: BED JSON entries with reference alleles.
    """
    bed_entries = []

    for chrom, (start, end) in studied_regions.items():
        try:
            # Extract the chromosome sequence using extract_chromosome_sequence
            chromosome_seq = extract_chromosome_sequence(chrom.lstrip("Chr"), build_version)

            # Extract the reference allele for this region
            ref_allele = chromosome_seq[start - 1:end]  # Convert 1-based to 0-based for slicing
        except ValueError as e:
            print(f"Error processing region {chrom}:{start}-{end}: {e}")
            ref_allele = ""

        # Append BED entry
        bed_entries.append({
            "CHROM": chrom,
            "START": start - 1,  # Convert to 0-based for BED format
            "END": end,  # Keep as 1-based
            "REF_ALLELE": ref_allele,
        })

    return bed_entries


def main(gvcf_file_path, build_version):
    """
    Main function to parse the GVCF file and generate BED JSON with reference alleles.

    Args:
        gvcf_file_path (str): Path to the GVCF file.
        build_version (int): Reference genome build version (37 or 38).
    """
    print(f"Parsing GVCF file: {gvcf_file_path} with build version {build_version}")
    studied_regions, blocks = parse_gvcf(gvcf_file_path, build_version)

    # Generate BED JSON for studied regions
    bed_json = generate_bed_json(studied_regions, build_version)

    # Output the BED JSON to a file
    bed_output_file = "bed_output_with_ref_and_metadata.json"
    with open(bed_output_file, "w") as bed_file:
        json.dump(bed_json, bed_file, indent=2)
    print(f"BED JSON with reference alleles saved to {bed_output_file}")

    # Output the detailed blocks to a file
    blocks_output_file = "non_variant_blocks.json"
    with open(blocks_output_file, "w") as blocks_file:
        json.dump(blocks, blocks_file, indent=2)
    print(f"Detailed blocks saved to {blocks_output_file}")


# Example usage
gvcf_file_path = "/Users/ethankakavetsis/Desktop/genomics-operations/utilities/NA18870.chr20.GRCh37.g.vcf"
build_version = 37  # Specify the reference genome version (37 or 38)

if __name__ == "__main__":
    main(gvcf_file_path, build_version)
