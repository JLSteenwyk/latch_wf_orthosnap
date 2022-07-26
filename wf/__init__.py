"""
identify subgroups of single-copy orthologous genes in multi-copy gene families
"""

from enum import Enum
import glob
import os
from pathlib import Path
import subprocess
from typing import Optional, List

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir, file_glob
from latch.types.utils import _is_valid_url

class InparalogToKeep(Enum):
    none = "none"
    shortest_seq_len = "shortest_seq_len"
    median_seq_len = "median_seq_len"
    longest_seq_len = "longest_seq_len"
    shortest_branch_len = "shortest_branch_len"
    median_branch_len = "median_branch_len"
    longest_branch_len = "longest_branch_len"

# def file_glob(
#     pattern: str,
#     remote_directory: str,
#     target_dir: Optional[Path] = None
# ) -> List[LatchFile]:
#     """Constructs a list of LatchFiles from a glob pattern.
#     Convenient utility for passing collections of files between tasks. See
#     [nextflow's channels](https://www.nextflow.io/docs/latest/channel.html) or
#     [snakemake's wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards).
#     for similar functionality in other orchestration tools.
#     The remote location of each constructed LatchFile will be consructed by
#     appending the file name returned by the pattern to the directory
#     represented by the `remote_directory`.
#     Args:
#         pattern: A glob pattern to match a set of files, eg. '\*.py'. Will
#             resolve paths with respect to the working directory of the caller.
#         remote_directory: A valid latch URL pointing to a directory, eg.
#             latch:///foo. This _must_ be a directory and not a file.
#         target_dir: An optional Path object to define an alternate working
#             directory for path resolution
#     Returns:
#         A list of instantiated LatchFile objects.
#     Intended Use: ::
#         @small_task
#         def task():
#             ...
#             return file_glob("*.fastq.gz", "latch:///fastqc_outputs")
#     """

#     if not _is_valid_url(remote_directory):
#         return []

#     if target_dir is None:
#         wd = Path.cwd()
#     else:
#         wd = target_dir
#     matched = sorted(wd.glob(pattern))

#     return [LatchFile(str(file), remote_directory + file.name) for file in matched]

@small_task
def orthosnap_task(
    multi_copy_gene_tree: LatchFile,
    multi_copy_fasta_file: LatchFile,
    output_dir: LatchDir,
    inparalog_to_keep: InparalogToKeep,
    rooted: bool = False, 
    snap_trees: bool = False,
    support: Optional[float] = 80.0,
    occupancy: Optional[int] = 5,
) -> LatchDir: #List[LatchFile]:

    if inparalog_to_keep.value == "none":
        inparalog_to_keep = InparalogToKeep.longest_seq_len

    if support == None:
        support = 0.0

    # OrthoSNAP handling
    _orthosnap_cmd = [
        "orthosnap",
        "-t",
        multi_copy_gene_tree.local_path,
        "-f",
        multi_copy_fasta_file.local_path,
        "-s",
        str(support),
        "-o",
        str(occupancy),
        "-ip",
        inparalog_to_keep.value,
    ]

    if rooted: # if input tree is rooted
        _orthosnap_cmd.append("-r")

    if snap_trees: # if SNAP-OG trees are to be outputted
        _orthosnap_cmd.append("-st")

    local_dir = "/root/orthosnap_output/" #local directory to put output files in

    _orthosnap_cmd.append("-op")
    _orthosnap_cmd.append(local_dir)

    subprocess.run(_orthosnap_cmd)

    # return file_glob(f"{local_dir}/*.orthosnap.*", output_dir.remote_path) #, Path(multi_copy_fasta_file))
    return LatchDir(local_dir, output_dir.remote_path)

@workflow
def orthosnap(
    multi_copy_gene_tree: LatchFile,
    multi_copy_fasta_file: LatchFile,
    output_dir: LatchDir,
    inparalog_to_keep: InparalogToKeep,
    rooted: bool = False, 
    snap_trees: bool = False, 
    support: Optional[float] = 80.0,
    occupancy: Optional[int] = 5,
) -> LatchDir: #List[LatchFile]:
    """
    OrthoSNAP
    ----
    # OrthoSNAP: a tree splitting and pruning algorithm for retrieving single-copy orthologs from gene family trees
    ## About
    OrthoSNAP is a command-line tool for increasing the size of molecular evolution datasets that can then be used for
    diverse studies including phylogenomics and genome-wide surveys of selection.

    <br />

    Molecular evolution studies, such as phylogenomic studies and genome-wide surveys of selection, often rely on gene
    families of single-copy orthologs (SC-OGs). Large gene families with multiple homologs in one or more species are
    often ignored because identifying and retrieving SC-OGs nested within them is challenging. To address this issue
    and increase the number of markers used in molecular evolution studies, we developed OrthoSNAP, a software that
    uses a phylogenetic framework to simultaneously split gene families into SC-OGs and prune species-specific
    inparalogs. We term SC-OGs identified by OrthoSNAP as SNAP-OGs because they are identified using a splitting and
    pruning procedure (analogous to snapping branches on a tree). OrthoSNAP is useful for increasing the number of
    markers used in molecular evolution data matrices, a critical step for robustly inferring and exploring the tree of
    life.

    <br />

    ## Arguments
    ### Required
    -t, --tree <newick tree file>
        - input tree file in newick format
        - taxa name and gene name should be separate by a "|" symbol
            For example, "gene_a" from "species_a" should appear in the
            tree as "species_a|gene_a", and "gene_b" from "species_a"
            should appear in the tree as "species_a|gene_b", and so 
            on and so forth.
    -f, --fasta <fasta file>
        - input sequence file in fasta format 
        - taxa name and gene name should be formatted the same as in
            the tree file. Thus, "gene_a" from "species_a" should appear
            in the tree as "species_a|gene_a" and so on and so forth.
    
    ### Optional arguments
    -s, --support  <support>
        - support threshold for bipartition collapsing
        - all bipartitions will values less than the specified value 
            will be collapsed. For example, if the support threshold
            value is 80, all bipartitions with 79 or less support
            will be collapsed.
        - default value is 80 and is set for ultrafast bootstrap
            approximations. If bipartitions support was evaluated 
            using standard bootstrap, a common threshold to use is 70.
    -o, --occupancy  <occupancy>
        - occupancy threshold for minimum number of tips in orthologous
            subgroup. 
        - default value is 50 percent of the total number of taxa
        - values are rounded to the nearest integer. For example,
            if there are 15 taxa, the occupancy threshold will be 8.
    -r, --rooted
        - boolean argument for whether the input phylogeny is already rooted
        - if used, the input phylogeny is assumed to be rooted; if not,
            the tree will be midpoint rooted
    -st, --snap_trees
        - boolean argument for whether newick files of SNAP-OGs should also
            be outputted
        - if used, newick files of SNAP-OGs will be outputteds
    -ip, --inparalog_to_keep <shortest_seq_len,
                                median_seq_len,
                                longest_seq_len,
                                shortest_branch_len,
                                median_branch_len,
                                longest_branch_len>                                
        - specify how to determine which species-specific inparalog should be kept
        - the species-specific inparalog can be kept based on sequence length
            (shortest/median/longest_seq_len) or branch length based on tip-to-root
            distances (shortest/median/longest_branch_len)
        - by default, the longest sequence is kept following the standard approach
            in transcriptomics

    <br />

    ## Citation

    If you found OrthoSNAP useful, please cite *OrthoSNAP: a tree splitting and pruning algorithm for retrieving
    single-copy orthologs from gene family trees*. Steenwyk et al. 2021, bioRxiv. doi:
    [10.1101/2021.10.30.466607v1](https://www.biorxiv.org/content/10.1101/2021.10.30.466607v1).

    __metadata__:
        display_name: OrthoSNAP
        author:
            name: Jacob L. Steenwyk
            email: jlsteenwyk@gmail.com
            github: https://github.com/JLSteenwyk
            twitter: JLSteenwyk
        repository: https://github.com/JLSteenwyk/orthosnap
        license:
            id: MIT


    Args:

        multi_copy_gene_tree:
            Newick tree file of a multi-copy gene family
            __metadata__:
                display_name: "Input multi-copy gene family tree"
                appearance:
					comment: "Input multi-copy gene family tree"

        multi_copy_fasta_file:
            FASTA file of a multi-copy gene family
            __metadata__:
                display_name: "Input multi-copy gene family FASTA"
                appearance:
					comment: "Input multi-copy gene family FASTA"

        output_dir:
            Output directory
			__metadata__:
				display_name: "Output directory"
				appearance:
					comment: "Output directory"

        inparalog_to_keep:
            Determine which species-specific inparalog to keep
            __metadata__:
                display_name: "Determine which species-specific inparalog to keep"
                appearance:
					comment: "- shortest_seq_len,
                    - median_seq_len,
                    - longest_seq_len,
                    - shortest_branch_len,
                    - median_branch_len,
                    - longest_branch_len}"

        support:
            bipartition support value for collapsing low supported branches
            __metadata__:
                display_name: "Bipartition support threshold"
                appearance:
					comment: "Bipartition support threshold"

        occupancy:
            minimum occupancy threshold for sequences in a SNAP-OG
            __metadata__:
                display_name: "Minimum occupancy threshold"
                appearance:
					comment: "Minimum occupancy threshold"

        rooted:
            Boolean for whether the input tree is already rooted 
            __metadata__:
                display_name: "Rooted input tree"

        snap_trees:
            Boolean for whether SNAP-OG phylogenies should also be outputted 
            __metadata__:
                display_name: "Output Newick files of SNAP-OGs"

    """

    return orthosnap_task(
        multi_copy_gene_tree=multi_copy_gene_tree,
        multi_copy_fasta_file=multi_copy_fasta_file,
        output_dir=output_dir,
        inparalog_to_keep=inparalog_to_keep,
        support=support,
        occupancy=occupancy,
        rooted=rooted,
        snap_trees=snap_trees,
    )
