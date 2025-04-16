from Bio import SeqIO
import pandas as pd
import os
from itertools import chain
import random
import re


def get_random_start(seed, contig_length, n):
    random.seed(seed)
    return [random.randint(0, contig_length - 1) for _ in range(n)]


def split_fasta(fa, outdir):
    record_ids = []
    name = os.path.basename(fa).split(".fasta")[0]
    fasta_name = os.path.join(outdir, name)

    for record in SeqIO.parse(fa, "fasta"):
        SeqIO.write(record, fasta_name + "_" + record.id + ".fna", "fasta")
        record_ids.append(
            {"contig": record.id, "fasta": name, "contig_length": len(record.seq)}
        )
    return record_ids


cov_df = pd.read_csv(snakemake.input.cov, sep="\t", dtype={"tax_id": int, "seed": int})

os.mkdir(snakemake.output.dir)
id2fa = []
for fa in snakemake.input.fa:
    id2fa.append(split_fasta(fa, snakemake.output.dir))
id2fa = list(chain.from_iterable(id2fa))
contig_df = pd.DataFrame.from_records(id2fa)
df = pd.merge(contig_df, cov_df, how="left", on="fasta")

cols = ["samplename", "fasta", "contig", "contig_length", "tax_id", "seed", "cov_sim"]

if snakemake.params.circular:
    cols += ["n", "random_start", "rotate"]

    if "rotate" not in df.columns:
        df.loc[:, "rotate"] = [snakemake.params.rotate] * len(df)

    # Apply pattern-based rules if auto-detect is enabled
    if hasattr(snakemake.params, 'auto_detect_circular') and snakemake.params.auto_detect_circular:
        # Pattern for circular contigs (plasmid, chromosome, or ptg/utg ending with 'c')
        circular_pattern = re.compile(r'(plasmid|chromosome|[up]tg\d+c)', re.IGNORECASE)
        # Pattern for explicitly linear contigs
        linear_pattern = re.compile(r'[up]tg\d+l', re.IGNORECASE)
        
        # Process all contigs
        for idx, row in df.iterrows():
            contig_name = row["contig"]
            
            # Force linear for anything with 'l' suffix pattern regardless of other settings
            if linear_pattern.search(contig_name):
                df.loc[idx, "rotate"] = 1
            # Apply circular pattern only if not already matched as linear
            elif not circular_pattern.search(contig_name):
                df.loc[idx, "rotate"] = 1
                
    df["random_start"] = df.apply(
        lambda row: get_random_start(row["seed"], row["contig_length"], row["rotate"]),
        axis=1,
    )
    df_expanded = df.explode("random_start").reset_index(drop=True)
    df_expanded["n"] = df_expanded.groupby(["samplename", "contig"]).cumcount() + 1
    df = df_expanded
    df["cov_sim"] = df["cov_sim"] / df["rotate"]

df[cols].to_csv(snakemake.output.tsv, sep="\t", index=False)
