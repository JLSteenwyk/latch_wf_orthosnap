"""
identify subgroups of single-copy orthologous genes in multi-copy gene families
"""

from enum import Enum
import subprocess
from typing import Optional

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir

class InparalogToKeep(Enum):
    shortest_seq_len = "shortest_seq_len"
    median_seq_len = "median_seq_len"
    longest_seq_len = "longest_seq_len"
    shortest_branch_len = "shortest_branch_len"
    median_branch_len = "median_branch_len"
    longest_branch_len = "longest_branch_len"

@small_task
def orthosnap_task(
    multi_copy_gene_tree: LatchFile,
    multi_copy_fasta_file: LatchFile,
    output_dir: LatchDir,
    support: Optional[float] = 80.0,
    occupancy: Optional[float] = None,
    rooted: Optional[bool] = False, 
    snap_trees: Optional[bool] = False, 
    inparalog_to_keep: Optional[InparalogToKeep] = InparalogToKeep.longest_seq_len,
    ) -> LatchDir:

    local_dir = "/root/orthosnap_output/" #local directory to put output files in

    # iqtree command
    _orthosnap_cmd = [
        "orthosnap",
        "-t",
        multi_copy_gene_tree.local_path,
        "-f",
        multi_copy_fasta_file.local_path,
        "-s",
        support,
        "-o",
        occupancy,
        "-r",
        rooted,
        "-st",
        snap_trees,
        "-ip",
        inparalog_to_keep.value,
    ]

    subprocess.run(_orthosnap_cmd)
    return LatchDir(local_dir, output_dir.remote_path) # this returns the directory in which all output files are stored

@workflow
def orthosnap(
    multi_copy_gene_tree: LatchFile,
    multi_copy_fasta_file: LatchFile,
    output_dir: LatchDir,
    support: Optional[float] = 80.0,
    occupancy: Optional[float] = None,
    rooted: Optional[bool] = False, 
    snap_trees: Optional[bool] = False, 
    inparalog_to_keep: Optional[InparalogToKeep] = InparalogToKeep.longest_seq_len,
    ) -> LatchDir:
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
        display_name: Identify subgroups of single-copy groups of genes using OrthoSNAP.
        author: Jacob L. Steenwyk
            name: Jacob L. Steenwyk
            email: jlsteenwyk@gmail.com
            github: https://github.com/JLSteenwyk
        repository: https://github.com/JLSteenwyk/orthosnap
        license:
            id: GNU-GPL


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
                appearance:
					comment: "- True
                    - False"

        snap_trees:
            Boolean for whether SNAP-OG phylogenies should also be outputted 
            __metadata__:
                display_name: "Rooted input tree"
                appearance:
					comment: "- True
                    - False"

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
                    - longest_branch_len"

        output_dir:
            Output directory
			__metadata__:
				display_name: "Output directory"
				appearance:
					comment: "Output directory"

    """

    return orthosnap_task(
        multi_copy_gene_tree=multi_copy_gene_tree,
        multi_copy_fasta_file=multi_copy_fasta_file,
        output_dir=output_dir,
        support=support,
        occupancy=occupancy,
        rooted=rooted,
        snap_trees=snap_trees,
        inparalog_to_keep=inparalog_to_keep,
    )
