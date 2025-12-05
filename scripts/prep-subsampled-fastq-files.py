import subprocess
import argparse
import glob
import os

"""
Script to automate subsampling FASTQ reads with rasusa at different total number of reads.
"""

def get_filenames(filename):
    """
    Basic function to get the basename of the input FASTQ files depending on the extension.

    Args:
        filename: Input filename to get the base of
    
    Returns: 
        Base filename without the extension
    """
    # strip off the directory name before working with name
    name = os.path.basename(filename)

    if name.endswith(".fastq.gz"):
        base = name.replace(".fastq.gz", "")
    else:
        base = name.replace(".fastq", "")
    
    if base.endswith("_R1"):
        base = base.replace("_R1", "")
    else:
        base = base.replace("_R2", "")
    
    return base

def run_rasusa(input_dir, num_reads_list, output_dir):
    """
    Automates `rasusa reads` command given a list of total reads to randomly subsample down to.

    Args:
        input_dir: Input directory where FASTQ files are located.
        num_reads_list: Input comma-separated list of total reads to subsample down to.
        output_dir: Output directory to put the subsampled output FASTQ files into.
    """
    # create output directory if doesn't already exist
    os.makedirs(output_dir, exist_ok = True)

    # parse the input list for holding input set of read depths to be processed
    num_reads = [int(x) for x in num_reads_list.split(",")]

    fastq_files = (
        glob.glob(os.path.join(input_dir, "*.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*.fastq"))
    )

    # group reads into pairs dictionary so isn't iterating through every single file
    pairs = {}
    for file in fastq_files:
        base = get_filenames(file)
        pairs.setdefault(base, {})
        if "_R1" in file:
            pairs[base]["R1"] = file
        elif "_R2" in file:
            pairs[base]["R2"] = file
    
    # iterate per sample
    for base, reads in pairs.items():
        r1 = reads.get("R1")
        r2 = reads.get("R2")

    # iterate per input read_depth and call the rasusa command
        for read_depth in num_reads:
            output_fastq_file_1 = os.path.join(output_dir, f"{base}_n{read_depth}_R1.fastq.gz")
            output_fastq_file_2 = os.path.join(output_dir, f"{base}_n{read_depth}_R2.fastq.gz")

            # rasusa command
            cmd = [
                'rasusa', 'reads',
                str(r1), str(r2),
                "--num", str(read_depth),
                "--seed", "1", # random seed
                "--output-type", "g", # output as .gzip
                "-o", str(output_fastq_file_1),
                "-o", str(output_fastq_file_2)
            ]
            
            print(f"Running Rasusa on sample {base} | Depth={read_depth}")
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            if os.path.exists(output_fastq_file_1) and os.path.exists(output_fastq_file_2):
                print(f"Wrote {base} sample subsampled at {read_depth} to {output_dir}/!")

def parse_arguments():
    parser = argparse.ArgumentParser(description = "Randomly subsample reads in a FASTQ file with Rasusua")
    parser.add_argument("--input_dir", help="Input directory containing FASTQ files ending in .fastq or .gz")
    parser.add_argument("--num_reads_list", type = str, help="Comma-separated list of read values to subset to")
    parser.add_argument("--output_dir", help="Output directory for subsampled FASTQ files")
    return parser.parse_args()

def main():
    args = parse_arguments()
    run_rasusa(args.input_dir, args.num_reads_list, args.output_dir)

if __name__ =="__main__":
    main()