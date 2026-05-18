#!/usr/bin/env python3
"""
PhyloArk: Taxonomy-aware downsampling

Selects a target number of genomes while preserving the descendant count
structure of an induced NCBI taxonomy tree. Within terminal bins, genomes
are ranked using a CheckM derived quality score.

Inputs:
  - metadata table with accession and taxid columns
  - NCBI taxdump directory containing nodes.dmp, names.dmp, merged.dmp

Outputs:
  - <out_prefix>.selected.tsv
  - <out_prefix>.selected_accessions.txt
  - <out_prefix>.clade_summary.tsv
  - <out_prefix>.skipped.tsv
  - <out_prefix>.params.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

NA_VALUES = {"", "NA", "N/A", "NaN", "nan", "null", "NULL", "None", "none"}
DIRECT_PREFIX = "__DIRECT_GENOMES_AT_NODE__"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select taxonomy-structure-preserving genome representatives."
    )
    parser.add_argument("--genomes", required=True, help="Genome metadata table, TSV by default.")
    parser.add_argument("--taxdump", required=True, help="Directory containing NCBI nodes.dmp and names.dmp.")
    parser.add_argument("--target", type=int, default=1000, help="Target number of genomes to select.")
    parser.add_argument(
        "--root-taxid",
        default=None,
        help="Root taxid to use as the top of the induced tree. Genomes outside this root are skipped. (OPTIONAL)",
    )
    parser.add_argument("--out-prefix", default="phyloark_stability_sample", help="Output file prefix.")
    parser.add_argument("--delimiter", default="\t", help="Input/output delimiter. Default: tab.")
    parser.add_argument("--accession-col", default="accession", help="Accession column name.")
    parser.add_argument("--taxid-col", default="taxid", help="Taxid column name.")
    parser.add_argument(
        "--score-col",
        default=None,
        help="Precomputed quality score column ranking genomes. Default score is computed from CheckM completeness and contamination if none is provided. (OPTIONAL)",
    )
    parser.add_argument("--complete-col", default="checkM_complete", help="CheckM completeness col name.")
    parser.add_argument("--contam-col", default="checkM_contam", help="CheckM contamination col name.")
    parser.add_argument(
        "--contam-weight",
        type=float,
        default=5.0,
        help="Penalty multiplier for contamination when computing CheckM quality score. Default: 5.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.5,
        help="Quota smoothing exponent. 1=proportional, 0=equal child clades, 0.5=sqrt smoothing. Default: 0.5.",
    )
    parser.add_argument(
        "--no-min-one",
        action="store_true",
        help="Do not force one genome per child bin when local budget allows it.",
    )
    parser.add_argument(
        "--rank-keep",
        default=None,
        help=(
            "Optional comma-separated rank whitelist for internal nodes, e.g. "
            "family,genus,species. The root and each genome's own taxid are always kept. "
            "Default: keep all NCBI taxonomy nodes in the induced tree."
        ),
    )
    parser.add_argument(
        "--include-lineage-names",
        action="store_true",
        help="Include lineage names in the selected output.",
    )
    return parser.parse_args()


def dmp_fields(line: str) -> List[str]:
    return [part.strip() for part in line.rstrip("\n").split("|")]


def read_nodes(nodes_path: Path) -> Tuple[Dict[int, int], Dict[int, str]]:
    parent: Dict[int, int] = {}
    rank: Dict[int, str] = {}
    with nodes_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = dmp_fields(line)
            try:
                taxid = int(fields[0])
                parent_taxid = int(fields[1])
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Could not parse nodes.dmp line: {line!r}") from exc
            parent[taxid] = parent_taxid
            rank[taxid] = fields[2] if len(fields) > 2 else ""
    return parent, rank


def read_names(names_path: Path) -> Dict[int, str]:
    names: Dict[int, str] = {}
    if not names_path.exists():
        return names
    with names_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = dmp_fields(line)
            if len(fields) < 4:
                continue
            try:
                taxid = int(fields[0])
            except ValueError:
                continue
            name_txt = fields[1]
            name_class = fields[3]
            if name_class == "scientific name":
                names[taxid] = name_txt
    return names


def read_merged(merged_path: Path) -> Dict[int, int]:
    merged: Dict[int, int] = {}
    if not merged_path.exists():
        return merged
    with merged_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = dmp_fields(line)
            if len(fields) < 2:
                continue
            try:
                old_taxid = int(fields[0])
                new_taxid = int(fields[1])
            except ValueError:
                continue
            merged[old_taxid] = new_taxid
    return merged


def parse_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text in NA_VALUES:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_taxid(value: object) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if text in NA_VALUES:
        return None
    try:
        return int(text)
    except ValueError:
        try:
            as_float = float(text)
            if as_float.is_integer():
                return int(as_float)
        except ValueError:
            return None
    return None


def read_table(path: Path, delimiter: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def resolve_taxid(taxid: int, parent: Dict[int, int], merged: Dict[int, int]) -> Optional[int]:
    seen = set()
    cur = taxid
    while cur in merged and cur not in seen:
        seen.add(cur)
        cur = merged[cur]
    if cur in parent:
        return cur
    return None


def path_to_root(taxid: int, parent: Dict[int, int]) -> Optional[List[int]]:
    path: List[int] = []
    seen = set()
    cur = taxid
    while True:
        if cur in seen:
            return None
        seen.add(cur)
        path.append(cur)
        if cur not in parent:
            return None
        par = parent[cur]
        if par == cur:
            return path
        cur = par


def lca_from_leaf_to_root_paths(paths: Sequence[List[int]]) -> int:
    if not paths:
        raise ValueError("Cannot compute LCA of zero paths.")
    root_to_leaf = [list(reversed(path)) for path in paths]
    min_len = min(len(path) for path in root_to_leaf)
    lca = root_to_leaf[0][0]
    for idx in range(min_len):
        candidate = root_to_leaf[0][idx]
        if all(path[idx] == candidate for path in root_to_leaf):
            lca = candidate
        else:
            break
    return lca


def lineage_root_to_leaf(
    leaf_to_root: List[int], root_taxid: int
) -> Optional[List[int]]:
    try:
        root_index = leaf_to_root.index(root_taxid)
    except ValueError:
        return None
    return list(reversed(leaf_to_root[: root_index + 1]))


def filter_lineage_by_rank(
    lineage: List[int], rank: Dict[int, str], keep_ranks: Optional[set]
) -> List[int]:
    if not keep_ranks:
        return lineage
    root = lineage[0]
    leaf = lineage[-1]
    out: List[int] = []
    for node in lineage:
        if node == root or node == leaf or rank.get(node, "") in keep_ranks:
            if not out or out[-1] != node:
                out.append(node)
    return out


def safe_score(row: Dict[str, str], args: argparse.Namespace) -> Tuple[float, float, float]:
    completeness = parse_float(row.get(args.complete_col))
    contamination = parse_float(row.get(args.contam_col))
    if contamination is None:
        contamination = 0.0

    if args.score_col:
        score = parse_float(row.get(args.score_col))
        if score is None:
            score = float("-inf")
    else:
        if completeness is None:
            score = float("-inf")
        else:
            score = completeness - args.contam_weight * contamination

    comp_for_sort = completeness if completeness is not None else float("-inf")
    contam_for_sort = contamination if contamination is not None else float("inf")
    return score, comp_for_sort, contam_for_sort


def make_added_headers(existing: Sequence[str], requested: Sequence[str]) -> List[str]:
    existing_set = set(existing)
    result = []
    for name in requested:
        candidate = name
        counter = 2
        while candidate in existing_set or candidate in result:
            candidate = f"{name}_{counter}"
            counter += 1
        result.append(candidate)
    return result


def allocate_budget(
    items: List[Tuple[object, int, int]],
    budget: int,
    alpha: float,
    min_one: bool,
) -> Dict[object, int]:
    """Allocate integer budget across capped items.

    Each item is (key, capacity, full_count). full_count controls the quota weight;
    capacity caps the number of selectable genomes. Returns exact allocation if
    budget <= total capacity.
    """
    clean_items = [(key, int(cap), int(count)) for key, cap, count in items if int(cap) > 0]
    total_capacity = sum(cap for _, cap, _ in clean_items)
    budget = min(int(budget), total_capacity)
    alloc: Dict[object, int] = {key: 0 for key, _, _ in clean_items}
    if budget <= 0 or not clean_items:
        return alloc

    remaining = budget
    if min_one and budget >= len(clean_items):
        for key, cap, _ in clean_items:
            if cap > 0:
                alloc[key] = 1
                remaining -= 1

    # Largest-remainder allocation with capacity caps
        available = [(key, cap, count) for key, cap, count in clean_items if alloc[key] < cap]
        if not available:
            break

        if alpha == 0:
            weights = {key: 1.0 for key, _, _ in available}
        else:
            weights = {key: max(float(count), 1.0) ** alpha for key, _, count in available}
        total_weight = sum(weights.values())
        if total_weight <= 0:
            weights = {key: 1.0 for key, _, _ in available}
            total_weight = float(len(available))

        raw = {key: remaining * weights[key] / total_weight for key, _, _ in available}
        floor_add: Dict[object, int] = {}
        remainders: Dict[object, float] = {}
        added = 0

        for key, cap, _ in available:
            room = cap - alloc[key]
            floor_value = int(math.floor(raw[key]))
            add = min(room, floor_value)
            floor_add[key] = add
            remainders[key] = raw[key] - floor_value
            added += add

        for key, add in floor_add.items():
            alloc[key] += add
        remaining -= added

        if remaining <= 0:
            break

        ranked = sorted(
            [(key, cap, count) for key, cap, count in available if alloc[key] < cap],
            key=lambda x: (remainders.get(x[0], 0.0), x[2], x[1], str(x[0])),
            reverse=True,
        )

        placed = 0
        for key, _, _ in ranked:
            if remaining <= 0:
                break
            alloc[key] += 1
            remaining -= 1
            placed += 1

        if placed == 0:
            break

    return alloc


def main() -> int:
    args = parse_args()
    if args.target <= 0:
        raise ValueError("--target must be positive.")
    if args.alpha < 0:
        raise ValueError("--alpha must be >= 0.")

    taxdump = Path(args.taxdump)
    nodes_path = taxdump / "nodes.dmp"
    names_path = taxdump / "names.dmp"
    merged_path = taxdump / "merged.dmp"
    if not nodes_path.exists():
        raise FileNotFoundError(f"Missing required file: {nodes_path}")

    parent, rank = read_nodes(nodes_path)
    names = read_names(names_path)
    merged = read_merged(merged_path)

    headers, rows = read_table(Path(args.genomes), args.delimiter)
    for required in [args.taxid_col, args.accession_col]:
        if required not in headers:
            raise ValueError(f"Missing required column: {required}")
    if args.score_col and args.score_col not in headers:
        raise ValueError(f"--score-col not found in table: {args.score_col}")
    if not args.score_col:
        for required in [args.complete_col, args.contam_col]:
            if required not in headers:
                raise ValueError(f"Missing required column for score calculation: {required}")

    keep_ranks = None
    if args.rank_keep:
        keep_ranks = {x.strip() for x in args.rank_keep.split(",") if x.strip()}

    root_taxid = parse_taxid(args.root_taxid) if args.root_taxid is not None else None
    if root_taxid is not None:
        resolved_root = resolve_taxid(root_taxid, parent, merged)
        if resolved_root is None:
            raise ValueError(f"Root taxid {root_taxid} is not present in nodes.dmp or merged.dmp.")
        root_taxid = resolved_root

    row_meta: Dict[int, Dict[str, object]] = {}
    valid_leaf_paths: Dict[int, List[int]] = {}
    skipped: List[Dict[str, object]] = []

    for row_id, row in enumerate(rows):
        raw_taxid = parse_taxid(row.get(args.taxid_col))
        accession = row.get(args.accession_col, "")
        if raw_taxid is None:
            skipped.append({"row_id": row_id, "accession": accession, "taxid": row.get(args.taxid_col, ""), "taxid_used": "", "reason": "missing_or_invalid_taxid"})
            continue
        used_taxid = resolve_taxid(raw_taxid, parent, merged)
        if used_taxid is None:
            skipped.append({"row_id": row_id, "accession": accession, "taxid": raw_taxid, "taxid_used": "", "reason": "taxid_not_found"})
            continue
        leaf_path = path_to_root(used_taxid, parent)
        if leaf_path is None:
            skipped.append({"row_id": row_id, "accession": accession, "taxid": raw_taxid, "taxid_used": used_taxid, "reason": "lineage_error"})
            continue
        valid_leaf_paths[row_id] = leaf_path
        score, completeness, contamination = safe_score(row, args)
        row_meta[row_id] = {
            "raw_taxid": raw_taxid,
            "taxid_used": used_taxid,
            "score": score,
            "completeness": completeness,
            "contamination": contamination,
            "accession": accession,
        }

    if not valid_leaf_paths:
        raise ValueError("No valid genomes remained after taxid parsing.")

    if root_taxid is None:
        root_taxid = lca_from_leaf_to_root_paths(list(valid_leaf_paths.values()))

    children: Dict[int, set] = defaultdict(set)
    genomes_at_node: Dict[int, List[int]] = defaultdict(list)
    lineage_by_row: Dict[int, List[int]] = {}
    all_nodes = set([root_taxid])

    for row_id, leaf_path in valid_leaf_paths.items():
        lineage = lineage_root_to_leaf(leaf_path, root_taxid)
        if lineage is None:
            meta = row_meta[row_id]
            skipped.append({
                "row_id": row_id,
                "accession": meta["accession"],
                "taxid": meta["raw_taxid"],
                "taxid_used": meta["taxid_used"],
                "reason": "outside_root_taxid",
            })
            continue
        lineage = filter_lineage_by_rank(lineage, rank, keep_ranks)
        lineage_by_row[row_id] = lineage
        all_nodes.update(lineage)
        for par, child in zip(lineage, lineage[1:]):
            children[par].add(child)
        genomes_at_node[lineage[-1]].append(row_id)

    if not lineage_by_row:
        raise ValueError("No genomes remained inside the selected root taxid.")

    # Convert child sets to sorted lists for deterministic traversal.
    sorted_children: Dict[int, List[int]] = {
        node: sorted(child_set, key=lambda x: (names.get(x, ""), x))
        for node, child_set in children.items()
    }

    sys.setrecursionlimit(max(10000, len(all_nodes) + 1000))
    desc_count_cache: Dict[int, int] = {}

    def desc_count(node: int) -> int:
        if node in desc_count_cache:
            return desc_count_cache[node]
        total = len(genomes_at_node.get(node, []))
        for child in sorted_children.get(node, []):
            total += desc_count(child)
        desc_count_cache[node] = total
        return total

    total_genomes = desc_count(root_taxid)
    target = min(args.target, total_genomes)
    min_one = not args.no_min_one

    def genome_sort_key(row_id: int):
        meta = row_meta[row_id]
        # Descending score and completeness, ascending contamination and accession.
        return (
            -float(meta["score"]),
            -float(meta["completeness"]),
            float(meta["contamination"]),
            str(meta["accession"]),
            row_id,
        )

    def select_from_node(node: int, budget: int) -> List[int]:
        budget = min(budget, desc_count(node))
        if budget <= 0:
            return []

        direct_rows = genomes_at_node.get(node, [])
        child_nodes = sorted_children.get(node, [])

        if not child_nodes:
            return sorted(direct_rows, key=genome_sort_key)[:budget]

        items: List[Tuple[object, int, int]] = []
        direct_key = (DIRECT_PREFIX, node)
        if direct_rows:
            items.append((direct_key, len(direct_rows), len(direct_rows)))
        for child in child_nodes:
            count = desc_count(child)
            items.append((child, count, count))

        alloc = allocate_budget(items, budget, args.alpha, min_one)
        selected: List[int] = []

        direct_budget = alloc.get(direct_key, 0)
        if direct_budget:
            selected.extend(sorted(direct_rows, key=genome_sort_key)[:direct_budget])

        for child in child_nodes:
            child_budget = alloc.get(child, 0)
            if child_budget:
                selected.extend(select_from_node(child, child_budget))

        # Defensive fill, should rarely be needed. Keeps exact target if rounding/caps
        # leave slack because of pathological input.
        if len(selected) < budget:
            already = set(selected)
            remaining_rows = []
            for row_id, lineage in lineage_by_row.items():
                if row_id not in already and node in lineage:
                    remaining_rows.append(row_id)
            selected.extend(sorted(remaining_rows, key=genome_sort_key)[: budget - len(selected)])
        return selected[:budget]

    selected_ids = select_from_node(root_taxid, target)
    selected_set = set(selected_ids)

    selected_count: Counter = Counter()
    full_count: Counter = Counter()
    for row_id, lineage in lineage_by_row.items():
        for node in lineage:
            full_count[node] += 1
        if row_id in selected_set:
            for node in lineage:
                selected_count[node] += 1

    # Depth map from induced tree.
    depth = {root_taxid: 0}
    stack = [root_taxid]
    while stack:
        node = stack.pop()
        for child in sorted_children.get(node, []):
            depth[child] = depth[node] + 1
            stack.append(child)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    added_requested = [
        "_phyloark_taxid_used",
        "_phyloark_checkm_score",
        "_phyloark_lineage_taxids",
    ]
    if args.include_lineage_names:
        added_requested.append("_phyloark_lineage_names")
    added_headers = make_added_headers(headers, added_requested)

    selected_path = Path(f"{out_prefix}.selected.tsv")
    with selected_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers + added_headers, delimiter=args.delimiter, extrasaction="ignore")
        writer.writeheader()
        for row_id in sorted(selected_ids, key=lambda rid: (lineage_by_row[rid], genome_sort_key(rid))):
            row = dict(rows[row_id])
            lineage = lineage_by_row[row_id]
            values = {
                added_headers[0]: row_meta[row_id]["taxid_used"],
                added_headers[1]: "-inf" if math.isinf(float(row_meta[row_id]["score"])) else f"{float(row_meta[row_id]['score']):.6f}",
                added_headers[2]: ";".join(str(x) for x in lineage),
            }
            if args.include_lineage_names:
                values[added_headers[3]] = ";".join(names.get(x, "") for x in lineage)
            row.update(values)
            writer.writerow(row)

    acc_path = Path(f"{out_prefix}.selected_accessions.txt")
    with acc_path.open("w") as handle:
        for row_id in selected_ids:
            handle.write(str(row_meta[row_id]["accession"]) + "\n")

    taxid_path = Path(f"{out_prefix}.selected_taxids.txt")
    with taxid_path.open("w") as handle:
        for row_id in selected_ids:
            handle.write(str(row_meta[row_id]["taxid_used"]) + "\n")

    summary_path = Path(f"{out_prefix}.clade_summary.tsv")
    summary_headers = [
        "taxid",
        "name",
        "rank",
        "depth",
        "direct_genomes",
        "children",
        "full_count",
        "selected_count",
        "full_prop",
        "selected_prop",
        "abs_prop_diff",
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_headers, delimiter="\t")
        writer.writeheader()
        for node in sorted(all_nodes, key=lambda x: (depth.get(x, 10**9), names.get(x, ""), x)):
            fc = full_count[node]
            sc = selected_count[node]
            full_prop = fc / total_genomes if total_genomes else 0.0
            sample_prop = sc / target if target else 0.0
            writer.writerow({
                "taxid": node,
                "name": names.get(node, ""),
                "rank": rank.get(node, ""),
                "depth": depth.get(node, ""),
                "direct_genomes": len(genomes_at_node.get(node, [])),
                "children": len(sorted_children.get(node, [])),
                "full_count": fc,
                "selected_count": sc,
                "full_prop": f"{full_prop:.8f}",
                "selected_prop": f"{sample_prop:.8f}",
                "abs_prop_diff": f"{abs(full_prop - sample_prop):.8f}",
            })

    skipped_path = Path(f"{out_prefix}.skipped.tsv")
    skipped_headers = ["row_id", "accession", "taxid", "taxid_used", "reason"]
    with skipped_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=skipped_headers, delimiter="\t")
        writer.writeheader()
        for item in skipped:
            writer.writerow(item)

    nonroot_nodes = [node for node in all_nodes if node != root_taxid]
    l1_distortion = sum(
        abs((full_count[node] / total_genomes) - (selected_count[node] / target))
        for node in nonroot_nodes
    ) if target else 0.0
    mean_abs_distortion = l1_distortion / len(nonroot_nodes) if nonroot_nodes else 0.0

    params = {
        "genomes": str(Path(args.genomes)),
        "taxdump": str(taxdump),
        "root_taxid": root_taxid,
        "root_name": names.get(root_taxid, ""),
        "target_requested": args.target,
        "target_used": target,
        "total_valid_genomes_under_root": total_genomes,
        "selected_genomes": len(selected_ids),
        "skipped_rows": len(skipped),
        "alpha": args.alpha,
        "min_one_per_child_when_possible": min_one,
        "rank_keep": sorted(keep_ranks) if keep_ranks else None,
        "score_col": args.score_col,
        "complete_col": args.complete_col,
        "contam_col": args.contam_col,
        "contam_weight": args.contam_weight,
        "l1_distortion_nonroot_nodes": l1_distortion,
        "mean_abs_distortion_nonroot_nodes": mean_abs_distortion,
        "outputs": {
            "selected_table": str(selected_path),
            "selected_accessions": str(acc_path),
            "selected_taxids": str(taxid_path),
            "clade_summary": str(summary_path),
            "skipped": str(skipped_path),
        },
    }
    params_path = Path(f"{out_prefix}.params.json")
    with params_path.open("w") as handle:
        json.dump(params, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Root: {root_taxid} {names.get(root_taxid, '')}".rstrip())
    print(f"Valid genomes under root: {total_genomes}")
    print(f"Selected genomes: {len(selected_ids)}")
    print(f"Skipped rows: {len(skipped)}")
    print(f"Mean abs distortion across non-root nodes: {mean_abs_distortion:.8f}")
    print(f"Selected table: {selected_path}")
    print(f"Clade summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(1)
