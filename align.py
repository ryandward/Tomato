import os
import argparse
import pysam
import subprocess
from collections import defaultdict
from rich.console import Console
from rich.highlighter import ReprHighlighter
from rich.progress import (
    Progress,
    TextColumn,
    BarColumn,
    TimeElapsedColumn,
    TransferSpeedColumn,
)
from multiprocessing import Process

console = Console()

# look https://github.com/nf-core/atacseq?tab=readme-ov-file for the guidelines


def run_bowtie2(fwd_read, rev_read, bowtie2_index, output_bam):
    count = 0
    bowtie2_command = [
        "bowtie2",
        "-p6",
        "--very-fast",
        "-x",
        bowtie2_index,
        "-1",
        fwd_read,
        "-2",
        rev_read,
    ]
    with subprocess.Popen(
        bowtie2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ) as bowtie2_process:
        with pysam.AlignmentFile(
            bowtie2_process.stdout, "rb", check_sq=False
        ) as samfile:
            with pysam.AlignmentFile(
                output_bam, "wb", header=samfile.header
            ) as bamfile:
                with Progress(
                    TextColumn(
                        "[progress.description]{task.description}",
                        highlighter=ReprHighlighter(),
                    ),
                    BarColumn(),
                    TextColumn(
                        "Reads: {task.fields[count]:,.0f}",
                        highlighter=ReprHighlighter(),
                    ),
                    TextColumn(
                        "(BAM size: {task.fields[output_size]:,.2f} MB)",
                        highlighter=ReprHighlighter(),
                    ),
                    TransferSpeedColumn(),
                ) as progress:
                    task = progress.add_task(
                        "Bowtie2...",
                        total=None,
                        output_size=0.0,
                        count=0,
                    )
                    for read in samfile:
                        if (
                            read.is_proper_pair
                            and not read.is_unmapped
                            and not read.mate_is_unmapped
                        ):
                            bamfile.write(read)
                            output_size = os.path.getsize(output_bam) / (
                                1024 * 1024
                            )  # Convert bytes to MB
                            count = count + 1
                            progress.update(task, output_size=output_size)
                            progress.update(task, count=count)
                            progress.update(task, completed=os.path.getsize(output_bam))
                    progress.update(task, total=count)


def sort_bam_process(bam_file, sorted_bam, threads=4):
    # Include the --threads option to utilize multiple cores for sorting
    pysam.sort("-o", sorted_bam, "-@", str(threads), bam_file)
    pysam.index(sorted_bam)


def sort_bam(bam_file, threads=4):
    sorted_bam = bam_file.replace(".bam", ".sorted.bam")

    process = Process(target=sort_bam_process, args=(bam_file, sorted_bam, threads))
    process.start()

    with Progress(
        TextColumn(
            "[progress.description]{task.description}", highlighter=ReprHighlighter()
        ),
        BarColumn(),
        TextColumn("Time elapsed:"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Sorting...", total=None)
        while process.is_alive():
            progress.refresh()  # Ensure the progress display is updated
        progress.update(task, total=1, advance=1)

    process.join()
    return sorted_bam


from pybedtools import BedTool


def parse_bam(bam_file, output_file):
    # Create a BedTool object from the BAM file
    bam = BedTool(bam_file)

    # Calculate the coverage and save it to a BEDGraph file
    bam.genome_coverage(bg=True, split=True, output=output_file)


def process_directory(directory, bowtie2_index):
    console.print(f"Processing directory: {os.path.abspath(directory)}")
    files = [f for f in os.listdir(directory) if f.endswith("_1.fastq")]
    experiments = set(f.split("_1.fastq")[0] for f in files)

    for exp in experiments:
        console.print(f"Processing experiment: [bold cyan]{exp}[/bold cyan]")
        fwd_read = os.path.join(directory, f"{exp}_1.fastq")
        rev_read = os.path.join(directory, f"{exp}_2.fastq")

        output_bam = os.path.join(directory, f"{exp}.bam")
        run_bowtie2(fwd_read, rev_read, bowtie2_index, output_bam)

        sorted_bam = sort_bam(output_bam)

        output_coords = os.path.join(directory, f"{exp}_coords.tsv")
        parse_bam(sorted_bam, output_coords)

        console.print(f"Finished processing experiment: [bold cyan]{exp}[/bold cyan]")


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
