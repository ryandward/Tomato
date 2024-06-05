import os
import argparse
import pysam
import subprocess
from collections import defaultdict
from rich.console import Console
from rich.highlighter import ReprHighlighter
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
)

console = Console()


def run_bowtie2(fwd_read, rev_read, bowtie2_index, output_bam):
    console.log("Starting bowtie2.")
    bowtie2_command = [
        "bowtie2",
        "-p10",
        "--very-fast",
        "-x",
        bowtie2_index,
        "-1",
        fwd_read,
        "-2",
        rev_read,
    ]
    samtools_view_command = ["samtools", "view", "-Sb", "-"]
    with subprocess.Popen(
        bowtie2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ) as bowtie2_process:
        with subprocess.Popen(
            samtools_view_command,
            stdin=bowtie2_process.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as samtools_view_process:
            with open(output_bam, "wb") as bam_file:
                with Progress(
                    SpinnerColumn(),
                    TextColumn(
                        "[progress.description]{task.description}",
                        highlighter=ReprHighlighter(),
                    ),
                    BarColumn(),
                    TextColumn(
                        "BAM file size: {task.fields[output_size]:,.2f} MB",
                        highlighter=ReprHighlighter(),
                    ),
                ) as progress:
                    task = progress.add_task(
                        "Aligning reads...",
                        total=None,
                        output_size=0.0,
                    )
                    while bowtie2_process.poll() is None:
                        bam_file.write(samtools_view_process.stdout.read(1024))
                        output_size = os.path.getsize(output_bam) / (
                            1024 * 1024
                        )  # Convert bytes to MB
                        progress.update(task, output_size=output_size)
                    bam_file.write(samtools_view_process.stdout.read())
                samtools_view_process.wait()
        bowtie2_process.wait()
    console.log("Bowtie2 alignment completed.")


def sort_bam(bam_file):
    console.log("Sorting BAM file.")
    sorted_bam = bam_file.replace(".bam", ".sorted.bam")
    subprocess.run(["samtools", "sort", "-o", sorted_bam, bam_file])
    subprocess.run(["samtools", "index", sorted_bam])
    return sorted_bam


def parse_and_shift_bam(bam_file, output_file):
    console.log("Parsing and shifting reads in BAM file.")
    bam = pysam.AlignmentFile(bam_file, "rb")
    coordinates_count = defaultdict(int)
    read_pairs_processed = 0

    total_reads = bam.mapped
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn(
            "{task.completed}/{task.total} pairs", highlighter=ReprHighlighter()
        ),
    ) as progress:
        task = progress.add_task("Processing reads...", total=round(total_reads / 2))

        for read in bam.fetch():
            if (
                read.is_proper_pair
                and not read.is_unmapped
                and not read.mate_is_unmapped
                and read.is_read1
            ):
                chrom = read.reference_name
                start = (
                    read.reference_start + 4
                    if not read.is_reverse
                    else read.reference_start - 5
                )
                end = (
                    read.next_reference_start + read.template_length - 5
                    if not read.mate_is_reverse
                    else read.next_reference_start + read.template_length + 4
                )

                for pos in range(start, end):
                    coordinates_count[(chrom, pos)] += 1

                read_pairs_processed += 1
                progress.update(task, advance=1)
                # if read_pairs_processed % 1000 == 0:
                #     console.print(
                #         f"{read_pairs_processed:,} read pairs have been processed.",
                #         end="\r",
                #     )

    with open(output_file, "w") as f:
        for (chrom, coord), count in coordinates_count.items():
            f.write(f"{chrom}\t{coord}\t{count}\n")


def process_directory(directory, bowtie2_index):
    console.log("Processing directory...")
    files = [f for f in os.listdir(directory) if f.endswith("_1.fastq")]
    experiments = set(f.split("_1.fastq")[0] for f in files)

    for exp in experiments:
        console.log(f"Processing experiment: {exp}")
        fwd_read = os.path.join(directory, f"{exp}_1.fastq")
        rev_read = os.path.join(directory, f"{exp}_2.fastq")

        output_bam = os.path.join(directory, f"{exp}.bam")
        run_bowtie2(fwd_read, rev_read, bowtie2_index, output_bam)

        sorted_bam = sort_bam(output_bam)

        output_coords = os.path.join(directory, f"{exp}_shifted_coords.tsv")
        parse_and_shift_bam(sorted_bam, output_coords)

        console.print(f"Finished processing experiment: {exp}")


def main():
    parser = argparse.ArgumentParser(
        description="Process paired-end reads, align them with Bowtie2, sort BAM files, and shift reads."
    )
    parser.add_argument(
        "directory", type=str, help="Directory containing the FASTQ files."
    )
    parser.add_argument(
        "bowtie2_index", type=str, help="Base name of the Bowtie2 index."
    )

    args = parser.parse_args()

    process_directory(args.directory, args.bowtie2_index)


if __name__ == "__main__":
    main()
