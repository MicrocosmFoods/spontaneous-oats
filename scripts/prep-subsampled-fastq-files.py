import subprocess
import argparse
import glob
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

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

def rasusa_single_sample_command(base, r1, r2, read_depth, output_dir):
    """
    Processes a single sample at a single read depth with `rasusa reads` to randomly subsample down an input set of paired-end FASTQ reads at given depth.

    Args:
        base: Base name of the sample
        r1: Path to R1 FASTQ file
        r2: Path to R2 FASTQ file
        read_depth: Number of reads to subsample to
        output_dir: Output directory to put the subsampled output FASTQ files into.
    """
    # define output files
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
    
    try:
        print(f"Running Rasusa on sample {base} | Depth={read_depth}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        if os.path.exists(output_fastq_file_1) and os.path.exists(output_fastq_file_2):
            print(f"Wrote {base} sample subsampled at {read_depth} to {output_dir}!")
            return(base, read_depth, True)
        else:
            return(base, read_depth, False)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {base} at depth {read_depth}: {e}")
        return(base, read_depth, False)

def run_rasusa(input_dir, num_reads_list, output_dir, num_threads =4):
    """
    Automates multi-processing of `rasusa reads` command given a list of total reads to randomly subsample down to.

    Args:
        input_dir: Input directory where FASTQ files are located
        num_reads_list: Input comma-separated list of total reads to susample down to
        output_dir: Output directory to write out subsampled output FASTQ files into
        num_threads: Number of parallel threads to use (default: 4)
    """
    # create output directory if doesn't exist already
    os.makedirs(output_dir, exist_ok = True)

    # parse input list of read depths
    num_reads = [int(x) for x in num_reads_list.split(",")]

    fastq_files = (
        glob.glob(os.path.join(input_dir, "*.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*.fastq"))
    )

    # group reads into pairs dictionary of dictionaries so isn't iterating through every single file and just a set of the pair
    pairs = {}
    for file in fastq_files:
        base = get_filenames(file)
        pairs.setdefault(base, {})
        if "_R1" in file:
            pairs[base]["R1"] = file
        elif "_R2" in file:
            pairs[base]["R2"] = file
    
    # create list of jobs depending on number of combo input fastq file pairs + read depth
    jobs = []
    for base, reads in pairs.items():
        r1 = reads.get("R1")
        r2 = reads.get("R2")
        for read_depth in num_reads:
            jobs.append((base, r1, r2, read_depth, output_dir))
    
    # process in parallel
    with ThreadPoolExecutor(max_workers = num_threads) as executor:
        futures = [executor.submit(rasusa_single_sample_command, *job) for job in jobs]

        # wait for jobs to complete
        for future in as_completed(futures):
            result = future.result() # output is base, read_depth, success status
 
def parse_arguments():
    parser = argparse.ArgumentParser(description = "Randomly subsample reads in a FASTQ file with Rasusua")
    parser.add_argument("--input_dir", help="Input directory containing FASTQ files ending in .fastq or .gz")
    parser.add_argument("--num_reads_list", type = str, help="Comma-separated list of read values to subset to")
    parser.add_argument("--output_dir", help="Output directory for subsampled FASTQ files")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads to use (default 4).")
    return parser.parse_args()

def main():
    args = parse_arguments()
    run_rasusa(args.input_dir, args.num_reads_list, args.output_dir, args.threads)

if __name__ =="__main__":
    main()