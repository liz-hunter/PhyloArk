#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
from collections import defaultdict


class DSU:
    def __init__(self):
        self.parent = []
        self.size = []

    def add(self):
        i = len(self.parent)
        self.parent.append(i)
        self.size.append(1)
        return i

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        return True


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_thresholds(s):
    vals = []
    for x in s.split(","):
        x = x.strip()
        if x:
            vals.append(float(x))
    return sorted(set(vals), reverse=True)


def is_header(fields):
    lowered = [x.lower() for x in fields]
    return "ani" in lowered and ("ref_file" in lowered or "query_file" in lowered)


def iter_skani_edges(path):
    """
    Yields:
      ref_file, query_file, ani, af_ref, af_query, ref_name, query_name
    """
    with open_text(path) as f:
        header = None
        col = None

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")

            if header is None:
                if is_header(fields):
                    header = fields
                    col = {name: i for i, name in enumerate(header)}
                    continue
                else:
                    header = []
                    col = {
                        "Ref_file": 0,
                        "Query_file": 1,
                        "ANI": 2,
                        "Align_fraction_ref": 3,
                        "Align_fraction_query": 4,
                        "Ref_name": 5,
                        "Query_name": 6,
                    }

            try:
                ref_file = fields[col.get("Ref_file", 0)]
                query_file = fields[col.get("Query_file", 1)]
                ani = float(fields[col.get("ANI", 2)])
                af_ref = float(fields[col.get("Align_fraction_ref", 3)])
                af_query = float(fields[col.get("Align_fraction_query", 4)])

                ref_name = ""
                query_name = ""
                if col.get("Ref_name", None) is not None and len(fields) > col["Ref_name"]:
                    ref_name = fields[col["Ref_name"]]
                if col.get("Query_name", None) is not None and len(fields) > col["Query_name"]:
                    query_name = fields[col["Query_name"]]

                yield ref_file, query_file, ani, af_ref, af_query, ref_name, query_name

            except Exception:
                print(f"WARNING: could not parse line, skipping: {line[:200]}", file=sys.stderr)
                continue


def load_genome_list(path, get_index):
    if path is None:
        return

    with open_text(path) as f:
        for line in f:
            genome = line.strip()
            if not genome or genome.startswith("#"):
                continue
            get_index(genome)


def component_summary(dsu, n):
    counts = defaultdict(int)
    for i in range(n):
        counts[dsu.find(i)] += 1

    sizes = list(counts.values())
    sizes.sort(reverse=True)

    n_clusters = len(sizes)
    largest = sizes[0] if sizes else 0
    singletons = sum(1 for x in sizes if x == 1)

    return n_clusters, largest, singletons


def iter_bits(mask):
    while mask:
        lsb = mask & -mask
        idx = lsb.bit_length() - 1
        yield idx
        mask ^= lsb


def looks_refseq(path):
    base = os.path.basename(path)
    return base.startswith("GCF_")


def write_outputs(
    outdir,
    threshold,
    min_af,
    genomes,
    name_by_file,
    dsu,
    degree,
    ani_sum,
):
    os.makedirs(outdir, exist_ok=True)

    clusters = defaultdict(list)
    for i in range(len(genomes)):
        clusters[dsu.find(i)].append(i)

    sorted_clusters = sorted(
        clusters.values(),
        key=lambda members: (-len(members), genomes[min(members)]),
    )

    reps_tsv = os.path.join(outdir, "representatives.tsv")
    reps_txt = os.path.join(outdir, "representatives.txt")
    clusters_tsv = os.path.join(outdir, "clusters.tsv")

    with open(reps_tsv, "w") as rep_out, open(reps_txt, "w") as rep_list, open(clusters_tsv, "w") as cl_out:
        rep_out.write(
            "cluster_id\trepresentative\tcluster_size\trep_degree_at_threshold\t"
            "rep_mean_ani_to_linked_members\tis_refseq_gcf\trep_name\n"
        )
        cl_out.write(
            "cluster_id\tgenome\trepresentative\tcluster_size\tgenome_name\n"
        )

        for cluster_id, members in enumerate(sorted_clusters, start=1):
            def rep_key(i):
                deg = degree[i]
                mean_ani = ani_sum[i] / deg if deg > 0 else 0.0
                name = name_by_file.get(genomes[i], "")
                is_complete = 1 if "complete genome" in name.lower() else 0
                is_gcf = 1 if looks_refseq(genomes[i]) else 0

                # Primary: most connected within cluster at threshold.
                # Tie-breakers: mean ANI, complete-genome wording, RefSeq GCF, stable file name.
                return (
                    deg,
                    mean_ani,
                    is_complete,
                    is_gcf,
                    os.path.basename(genomes[i]),
                )

            rep_i = max(members, key=rep_key)
            rep = genomes[rep_i]
            rep_name = name_by_file.get(rep, "")
            rep_deg = degree[rep_i]
            rep_mean_ani = ani_sum[rep_i] / rep_deg if rep_deg > 0 else 0.0
            is_gcf = 1 if looks_refseq(rep) else 0

            rep_list.write(f"{rep}\n")
            rep_out.write(
                f"{cluster_id}\t{rep}\t{len(members)}\t{rep_deg}\t"
                f"{rep_mean_ani:.5f}\t{is_gcf}\t{rep_name}\n"
            )

            for i in sorted(members, key=lambda x: genomes[x]):
                genome = genomes[i]
                genome_name = name_by_file.get(genome, "")
                cl_out.write(
                    f"{cluster_id}\t{genome}\t{rep}\t{len(members)}\t{genome_name}\n"
                )

    print(f"Wrote: {reps_txt}", file=sys.stderr)
    print(f"Wrote: {reps_tsv}", file=sys.stderr)
    print(f"Wrote: {clusters_tsv}", file=sys.stderr)
    print(f"Chosen ANI threshold: {threshold}", file=sys.stderr)
    print(f"Minimum aligned fraction: {min_af}", file=sys.stderr)
    print(f"Number of representative genomes: {len(sorted_clusters)}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Cluster skani sparse TSV output and choose representative genomes."
    )
    parser.add_argument(
        "--skani",
        required=True,
        help="skani triangle --sparse TSV output. May be .gz",
    )
    parser.add_argument(
        "--genomes",
        required=True,
        help="Original genome_paths.txt, one genome path per line.",
    )
    parser.add_argument(
        "--outdir",
        default="skani_reps",
        help="Output directory.",
    )
    parser.add_argument(
        "--thresholds",
        default="98,98.5,99,99.25,99.5,99.75",
        help="Comma-separated ANI thresholds to test if --ani is not supplied.",
    )
    parser.add_argument(
        "--ani",
        type=float,
        default=None,
        help="Use this ANI threshold directly instead of choosing by --target.",
    )
    parser.add_argument(
        "--min-af",
        type=float,
        default=80.0,
        help="Minimum aligned fraction required in BOTH directions.",
    )
    parser.add_argument(
        "--target",
        type=int,
        default=1000,
        help="Target number of clusters/representatives when sweeping thresholds.",
    )

    args = parser.parse_args()

    thresholds = [args.ani] if args.ani is not None else parse_thresholds(args.thresholds)

    if not thresholds:
        raise SystemExit("ERROR: no thresholds provided.")

    os.makedirs(args.outdir, exist_ok=True)

    genomes = []
    genome_to_idx = {}
    name_by_file = {}

    dsus = {t: DSU() for t in thresholds}

    def get_index(genome):
        if genome not in genome_to_idx:
            genome_to_idx[genome] = len(genomes)
            genomes.append(genome)
            for d in dsus.values():
                d.add()
        return genome_to_idx[genome]

    print("Loading genome list...", file=sys.stderr)
    load_genome_list(args.genomes, get_index)
    print(f"Loaded {len(genomes)} genomes from list.", file=sys.stderr)

    print("Streaming skani TSV and building clusters...", file=sys.stderr)
    n_edges_used = {t: 0 for t in thresholds}
    n_lines = 0

    for ref_file, query_file, ani, af_ref, af_query, ref_name, query_name in iter_skani_edges(args.skani):
        n_lines += 1

        if ref_file == query_file:
            continue

        min_af = min(af_ref, af_query)
        if min_af < args.min_af:
            continue

        i = get_index(ref_file)
        j = get_index(query_file)

        if ref_name:
            name_by_file.setdefault(ref_file, ref_name)
        if query_name:
            name_by_file.setdefault(query_file, query_name)

        for t in thresholds:
            if ani >= t:
                dsus[t].union(i, j)
                n_edges_used[t] += 1

        if n_lines % 5_000_000 == 0:
            print(f"Processed {n_lines:,} skani rows...", file=sys.stderr)

    print(f"Processed {n_lines:,} skani rows total.", file=sys.stderr)
    print(f"Total genomes seen: {len(genomes)}", file=sys.stderr)

    summary_path = os.path.join(args.outdir, "threshold_summary.tsv")
    rows = []

    with open(summary_path, "w") as out:
        out.write("ani_threshold\tmin_af\tclusters\tlargest_cluster\tsingletons\tedges_used\n")
        for t in thresholds:
            n_clusters, largest, singletons = component_summary(dsus[t], len(genomes))
            rows.append((t, n_clusters, largest, singletons, n_edges_used[t]))
            out.write(
                f"{t}\t{args.min_af}\t{n_clusters}\t{largest}\t{singletons}\t{n_edges_used[t]}\n"
            )

    print(f"Wrote: {summary_path}", file=sys.stderr)

    if args.ani is not None:
        chosen_t = args.ani
    else:
        chosen_t = min(rows, key=lambda x: (abs(x[1] - args.target), x[0]))[0]

    print(f"Chosen threshold: {chosen_t}", file=sys.stderr)

    chosen_dsu = dsus[chosen_t]

    # Second pass: score representative candidates by within-cluster degree and ANI sum.
    degree = [0] * len(genomes)
    ani_sum = [0.0] * len(genomes)

    print("Second pass: scoring candidate representatives...", file=sys.stderr)

    for ref_file, query_file, ani, af_ref, af_query, ref_name, query_name in iter_skani_edges(args.skani):
        if ref_file == query_file:
            continue

        if ani < chosen_t:
            continue

        min_af = min(af_ref, af_query)
        if min_af < args.min_af:
            continue

        if ref_file not in genome_to_idx or query_file not in genome_to_idx:
            continue

        i = genome_to_idx[ref_file]
        j = genome_to_idx[query_file]

        if chosen_dsu.find(i) != chosen_dsu.find(j):
            continue

        degree[i] += 1
        degree[j] += 1
        ani_sum[i] += ani
        ani_sum[j] += ani

    write_outputs(
        outdir=args.outdir,
        threshold=chosen_t,
        min_af=args.min_af,
        genomes=genomes,
        name_by_file=name_by_file,
        dsu=chosen_dsu,
        degree=degree,
        ani_sum=ani_sum,
    )


if __name__ == "__main__":
    main()