#!/usr/bin/env python3
"""
Author:  Stefan Neuenschwander
Purpose: Correct RBK assembly by RPB
Version: 0.1
"""

import argparse
import polars as pl


# -------------------------------------------------
# FUNCTIONS
# -------------------------------------------------

def get_args():
    """ Get the command-line arguments """

    parser = argparse.ArgumentParser(description='Correct RBK assembly by RPB')
    parser.add_argument('-r', '--pysam_rbk', help='Path to pysamstats RBK', type=str,
                        default='data/rbk_rbk.variation_13k.txt')
    parser.add_argument('-p', '--pysam_rpb', help='Path to pysamstats RPB', type=str,
                        default='data/rbk_rpb.variation_13k.txt')
    parser.add_argument('-o', '--outfile', help='Path to corrected assembly', type=str, default='out/corrected_fna.fna')
    return parser.parse_args()


def prep_pysam(pysam_path: str) -> pl.DataFrame:
    """ Read pysam and the specified columns """

    df = pl.read_csv(pysam_path, separator='\t')
    df = df.select([
        'chrom', 'pos', 'ref', 'reads_all', 'reads_fwd', 'reads_rev',
        'A', 'A_fwd', 'A_rev', 'T', 'T_fwd', 'T_rev',
        'G', 'G_fwd', 'G_rev', 'C', 'C_fwd', 'C_rev',
        'N', 'N_fwd', 'N_rev',
    ])
    return df


def add_realcov(comb: pl.DataFrame) -> pl.DataFrame:
    """ Add real coverage per position (do not count gaps as coverage) """

    comb = comb.with_columns(
        (pl.col("A_fwd") + pl.col("T_fwd") + pl.col("G_fwd") + pl.col("C_fwd") + pl.col("N_fwd")).alias(
            "reads_fwd_position"),
        (pl.col("A_rev") + pl.col("T_rev") + pl.col("G_rev") + pl.col("C_rev") + pl.col("N_rev")).alias(
            "reads_rev_position"),
        (pl.col("A") + pl.col("T") + pl.col("G") + pl.col("C") + pl.col("N")).alias("reads_all_position"),
    )
    return comb


def add_sb(comb: pl.DataFrame) -> pl.DataFrame:
    """ Add strand bias flag and score """

    # Get strand bias flag and score for each nucleotide
    sb_flag_A, sb_score_A = ratio_expr_with_score("A")
    sb_flag_T, sb_score_T = ratio_expr_with_score("T")
    sb_flag_G, sb_score_G = ratio_expr_with_score("G")
    sb_flag_C, sb_score_C = ratio_expr_with_score("C")

    # Add strand_bias (0/1) column based on majority base
    comb = comb.with_columns([
        pl.when(pl.col('majority_consensus') == 'A').then(sb_flag_A)
        .when(pl.col('majority_consensus') == 'T').then(sb_flag_T)
        .when(pl.col('majority_consensus') == 'G').then(sb_flag_G)
        .when(pl.col('majority_consensus') == 'C').then(sb_flag_C)
        .otherwise(None)
        .cast(pl.Int8)
        .alias('strand_bias'),

        # Add strand bias score column
        pl.when(pl.col('majority_consensus') == 'A').then(sb_score_A)
        .when(pl.col('majority_consensus') == 'T').then(sb_score_T)
        .when(pl.col('majority_consensus') == 'G').then(sb_score_G)
        .when(pl.col('majority_consensus') == 'C').then(sb_score_C)
        .otherwise(None)
        .alias('sb_score')
    ])

    return comb


def ratio_expr_with_score(nucleotide):
    """ Calculate SB ratio """

    eps = 1e-10

    fwd_ratio = pl.col(f"{nucleotide}_fwd") / (pl.col("reads_fwd_position") - pl.col(f"{nucleotide}_fwd") + eps)
    rev_ratio = pl.col(f"{nucleotide}_rev") / (pl.col("reads_rev_position") - pl.col(f"{nucleotide}_rev") + eps)

    fwd_expr = (fwd_ratio >= 2) & (rev_ratio <= 0.5)
    rev_expr = (fwd_ratio <= 0.5) & (rev_ratio >= 2)

    sb_score = pl.when(fwd_expr).then(fwd_ratio / (rev_ratio + eps)) \
        .when(rev_expr).then(rev_ratio / (fwd_ratio + eps)) \
        .otherwise(pl.lit(1.0))  # 1.0 = no strong bias

    sb_flag = fwd_expr ^ rev_expr  # boolean flag for SB

    return sb_flag, sb_score


def add_lowcov_flag(comb: pl.DataFrame) -> pl.DataFrame:
    """ Check for low coverage """

    comb = comb.with_columns(
        (pl.col("reads_all_position") <= 20).cast(pl.Int8).alias("low_cov")
    )
    return comb


def add_mixed_base_threshold(comb: pl.DataFrame) -> pl.DataFrame:
    """ Calculate the threshold (80% of total coverage) """

    comb = comb.with_columns(
        (pl.col("reads_all_position") * 0.8).round(1).alias("threshold")
    )
    return comb


def add_mixed_mayority_consensus_flag(comb: pl.DataFrame) -> pl.DataFrame:
    """ Check if the most common base passes the mixed bases threshold, add bol """

    comb = comb.with_columns(
        pl.when(pl.col("majority_consensus") == "A").then(pl.col("A"))
        .when(pl.col("majority_consensus") == "T").then(pl.col("T"))
        .when(pl.col("majority_consensus") == "G").then(pl.col("G"))
        .when(pl.col("majority_consensus") == "C").then(pl.col("C"))
        .when(pl.col("majority_consensus") == "N").then(pl.col("N"))
        .otherwise(None).alias("matches_to_major")
    )
    # Check if the value in the column specified by majority_consensus is smaller than 80% of the total coverage
    comb = comb.with_columns(
        (pl.col("matches_to_major") <= pl.col("threshold")).cast(pl.Int8).alias("ambiguous20plus")
    )
    return comb


def add_problem_count(comb: pl.DataFrame) -> pl.DataFrame:
    """ Count problems """

    comb = comb.with_columns(
        (pl.col("low_cov") + pl.col("strand_bias") + pl.col("ambiguous20plus")).alias("problems")
    )
    return comb


def join_and_correct(df_rbk: pl.DataFrame, df_rpb: pl.DataFrame) -> pl.DataFrame:
    """ Join RBK and RPB data, create corrected consensus """

    df_both = df_rbk.join(
        df_rpb,
        left_on=["chrom", "pos"],
        right_on=["chrom", "pos"],
        how="left",
        coalesce=True
    )

    # Label mismatches between original consensus and RPB majority consensus
    df_both = df_both.with_columns(
        (pl.col("ref") != pl.col("majority_consensus_right")).cast(pl.Int8).alias("mismatch")
    )

    # Correct consensus
    df_both = df_both.with_columns(
        pl.when((pl.col("mismatch") == 1) & (pl.col("problems") != 0) & (pl.col("problems_right") == 0)).then(
            pl.col("majority_consensus_right"))
        .when((pl.col("mismatch") == 1) & (pl.col("problems") != 0) & (pl.col("problems_right") != 0)).then(
            pl.lit("N"))
        .otherwise(pl.col("ref")).alias("final_consensus")
    )


    # Add info
    df_both = df_both.with_columns(
        pl.when((pl.col("mismatch") == 1) & (pl.col("strand_bias") == 1) & (pl.col("problems_right") == 0)).then(
            pl.lit('sbmm_correctable'))
        .when((pl.col("mismatch") == 1) & (pl.col("strand_bias") == 1) & (pl.col("problems_right") >= 1)).then(
            pl.lit('sbmm_not_correctable'))
        .when((pl.col("mismatch") == 0) & (pl.col("strand_bias") == 1) & (pl.col("problems_right") == 0)).then(
            pl.lit('sb_correctable'))
        .when((pl.col("mismatch") == 0) & (pl.col("strand_bias") == 1) & (pl.col("problems_right") >= 1)).then(
            pl.lit('sb_not_correctable'))
        .when((pl.col("mismatch") == 1) & (pl.col("strand_bias") == 0) & (pl.col("problems_right") == 0)).then(
            pl.lit('mm_correctable'))
        .when((pl.col("mismatch") == 1) & (pl.col("strand_bias") == 0) & (pl.col("problems_right") >= 1)).then(
            pl.lit('mm_not_correctable'))
        .otherwise(None).alias("comment")
    )
    return df_both


def pldf_to_fasta(df_both: pl.DataFrame, outfile: str) -> None:
    """ Group by chromosome and create fasta from majority_consensus """
    grouped_df = df_both.group_by("chrom").agg(
        pl.col("final_consensus").map_elements("".join, return_dtype=pl.Utf8).alias("sequence"))

    def create_multifasta(df):
        """ Function to create a multifasta string """
        multifasta = ""
        for row in df.iter_rows(named=True):
            header = f">{row['chrom']}\n"
            sequence = f"{row['sequence']}\n"
            multifasta += header + sequence
        return multifasta

    # Create the multifasta string
    multifasta_str = create_multifasta(grouped_df)

    # Write the multifasta string to a file
    with open(outfile, "w") as f:
        f.write(multifasta_str)


# -------------------------------------------------
# MAIN
# -------------------------------------------------

args = get_args()

### RBK
input_rbk = prep_pysam(args.pysam_rbk)

# Add majority consensus
input_rbk = input_rbk.with_columns(
    pl.col('ref').alias("majority_consensus")
)

# # Filter the DataFrame
# input_rbk_chr1 = input_rbk.filter(
# 	(pl.col('chrom') == 'contig_1') & 
# 	(pl.col('pos') >= 1)
# )

# process
processed_rbk = add_realcov(input_rbk)
processed_rbk = add_sb(processed_rbk)
processed_rbk = add_lowcov_flag(processed_rbk)
processed_rbk = add_mixed_base_threshold(processed_rbk)
processed_rbk = add_mixed_mayority_consensus_flag(processed_rbk)
processed_rbk = add_problem_count(processed_rbk)

# print RBK SB positions
processed_rbk_sb = processed_rbk.filter(
    (pl.col('strand_bias') == 1) &
    (pl.col('low_cov') == 0)
)
processed_rbk_sb_name = args.outfile.replace(".fna", "") + "_rbk_sb.csv"
# processed_rbk_sb.write_csv(processed_rbk_sb_name)
del [[input_rbk]]

### RPB
input_rpb = prep_pysam(args.pysam_rpb)

# Add majority consensus
input_rpb = input_rpb.with_columns(
    pl.struct(["A", "T", "G", "C", "N"]).map_elements(
        lambda x: None if all(v == 0 for v in x.values()) else max(x, key=x.get),
        return_dtype=pl.Utf8
    ).alias("majority_consensus")
)

# # Filter the DataFrame
# input_rpb_chr1 = input_rpb.filter(
# 	(pl.col('chrom') == 'contig_1') & 
# 	(pl.col('pos') >= 1)
# )

processed_rpb = add_realcov(input_rpb)
processed_rpb = add_sb(processed_rpb)
processed_rpb = add_lowcov_flag(processed_rpb)
processed_rpb = add_mixed_base_threshold(processed_rpb)
processed_rpb = add_mixed_mayority_consensus_flag(processed_rpb)
processed_rpb = add_problem_count(processed_rpb)
processed_rpb.write_csv('processed_rpb.csv')

del [[input_rpb]]

### RPB & RBK
df_both = join_and_correct(processed_rbk, processed_rpb)
pldf_to_fasta(df_both, args.outfile)

df_mm = df_both.filter(pl.col("mismatch") == 1)
df_mm.write_csv(args.outfile + '_mm_colls')
