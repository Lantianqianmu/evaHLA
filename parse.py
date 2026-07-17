#!/usr/bin/env python3
"""Merge HLA-HD and T1K genotyping results into a single flagged-style TSV.

Ports the logic from Code/evaHLA.R:
  1. Parse HLA-HD `*_HLA_HD_final.result.txt` files.
  2. Parse T1K `*_allele.tsv` files (dropping any allele with score < 60).
  3. For the 8 "simple" genes (A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1),
     prefer HLA-HD's call for a slot when it is "usable" (3-field, or a
     2-field call corroborated by / unopposed in T1K), otherwise fall back
     to T1K, picking whichever T1K allele is most different from the
     already-known slot when only one hla_hd slot needs filling.
  4. DRB3/4/5 is inferred entirely from T1K's raw DRB3/4/5 candidates,
     disambiguated using the DRB1 linkage rule (same table as HLA-HD's own
     resolution), plus an `unexpected` column listing any full-confidence
     T1K DRB3/4/5 candidate not selected by that rule.
  5. Homozygous pairs (identical alleles in both slots) are collapsed to
     report the allele once (slot 2 blanked).

Samples are merged purely by sample name (no `meta` map).

Usage:
    parse.py --hlahd DIR_OR_FILES --t1k DIR_OR_FILES --out OUT.tsv
    parse.py <sample> <hlahd_result.txt> <t1k_allele1.tsv> [<t1k_allele2.tsv> ...] --out OUT.tsv
"""
# python parse.py --sample H16 --hlahd /home/zeemeeuw/YangLab/TCR-seq/HLA/HLA_HD/H16_HLA_HD_final.result.txt --t1k /home/zeemeeuw/YangLab/TCR-seq/HLA/HLA_t1k/T1K_H16_R1_allele.tsv --out ./H16_parsed.tsv
# python parse.py --sample ${i} --hlahd /home/zeemeeuw/YangLab/TCR-seq/HLA/HLA_HD/${i}_HLA_HD_final.result.txt --t1k /home/zeemeeuw/YangLab/TCR-seq/HLA/HLA_t1k/T1K_${i}_R1_allele.tsv --out ./${i}_parsed.tsv
# python parse.py --sample AC12 --hlahd /home/zeemeeuw/YangLab/TCR-seq/HLA/HLA_HD/AC12_HLA_HD_final.result.txt --t1k /home/zeemeeuw/YangLab/TCR-seq/HLA/HLA_t1k/AC12_t1k_allele.tsv --out ./AC12_parsed.tsv

import argparse
import csv
import glob
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

SIMPLE_GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
OTHER_GENE_COLS = [f"{g}_{i}" for g in SIMPLE_GENES for i in (1, 2)]
DRB345_COLS = ["DRB3/4/5_1", "DRB3/4/5_2"]
OUTPUT_COLS = ["sample"] + OTHER_GENE_COLS + DRB345_COLS + ["unexpected", "QC_status", "QC_note"]

# DRB1 2-digit group -> expected DRB3/4/5 gene ("null" = no DRB3/4/5 partner)
DRB1_LINKAGE = {
    "01": "null", "03": "DRB3", "04": "DRB4", "07": "DRB4",
    "08": "null", "09": "DRB4", "10": "null", "11": "DRB3",
    "12": "DRB3", "13": "DRB3", "14": "DRB3", "15": "DRB5", "16": "DRB5",
}


# --------------------------------------------------------------------------
# Shared allele-string helpers
# --------------------------------------------------------------------------

def clean_allele(x: Optional[str]) -> Optional[str]:
    if x is None or x in ("-", "Not typed"):
        return None
    return re.sub(r"^HLA-", "", x)


def to_3field(x: Optional[str]) -> Optional[str]:
    """Truncate an allele string to its first 3 colon-separated fields."""
    if x is None:
        return None
    parts = x.split(":")
    return ":".join(parts[:3])


def field_count(x: Optional[str]) -> Optional[int]:
    if x is None:
        return None
    return len(x.split(":"))


def get_grp2(x: Optional[str]) -> Optional[str]:
    """First 2-digit group after '*' in an allele string, e.g. A*02:10 -> '02'."""
    if x is None:
        return None
    m = re.match(r"^[A-Za-z0-9]+\*(\d{2})", x)
    return m.group(1) if m else None


def diff_key(ref: str, cand: str) -> Tuple[int, int]:
    """(first differing field index, count of matching leading fields)
    between two colon-separated allele strings."""
    rf = ref.split(":")
    rc = cand.split(":")
    n = min(len(rf), len(rc))
    first_diff = n
    matches = 0
    for i in range(n):
        if rf[i] == rc[i]:
            matches += 1
        elif first_diff == n:
            first_diff = i
    return first_diff, matches


def pick_most_different(known: str, a1: Optional[str], a2: Optional[str],
                         default_a2: bool = True) -> Optional[str]:
    """Pick whichever of a1/a2 is most different from `known` (fewest
    shared leading fields); ties fall back to the positional default."""
    if a1 is None and a2 is None:
        return None
    if a1 is None:
        return a2
    if a2 is None:
        return a1
    k1 = diff_key(known, a1)
    k2 = diff_key(known, a2)
    if k1[0] != k2[0]:
        return a1 if k1[0] < k2[0] else a2
    if k1[1] != k2[1]:
        return a1 if k1[1] < k2[1] else a2
    return a2 if default_a2 else a1


def is_usable_hd(hd: Optional[str], own_data_raw: Optional[str],
                  other_data_raw: Optional[str]) -> bool:
    """An HLA-HD slot is usable if it's 3-field, or a 2-field call that
    either matches either T1K allele of the gene pair (corroborated) or
    has no competing T1K value at that slot (unopposed)."""
    hd_fc = field_count(hd)
    if hd_fc is None:
        return False
    if hd_fc == 3:
        return True
    if hd_fc != 2:
        return False
    own_fc = field_count(own_data_raw)
    other_fc = field_count(other_data_raw)
    corroborated = (own_fc == 2 and hd == own_data_raw) or \
                   (other_fc == 2 and hd == other_data_raw)
    unopposed = own_fc is None
    return corroborated or unopposed


def merge_drb345(gene_alleles: Dict[str, List[str]],
                  exp_genes: Tuple[Optional[str], Optional[str]]) -> Tuple[Optional[str], Optional[str]]:
    """Combine DRB3/DRB4/DRB5 candidate alleles into the final 2 observed
    slots. When more than 2 candidates are typed, use the DRB1-implied
    expected genes to pick the 2 most likely alleles; otherwise keep all
    typed alleles as-is."""
    entries = []
    for gene in ("DRB3", "DRB4", "DRB5"):
        for allele in gene_alleles.get(gene, []):
            entries.append((gene, allele))

    out = [None, None]
    if len(entries) <= 2:
        for i, (_, allele) in enumerate(entries):
            out[i] = allele
        return tuple(out)

    used = [False] * len(entries)
    slot = 0
    exp_clean = [g for g in exp_genes if g and g != "null"]
    for exp_gene in exp_clean:
        for i, (gene, allele) in enumerate(entries):
            if gene == exp_gene and not used[i]:
                out[slot] = allele
                used[i] = True
                slot += 1
                break
        if slot > 1:
            break

    for i, (_, allele) in enumerate(entries):
        if slot > 1:
            break
        if not used[i]:
            out[slot] = allele
            used[i] = True
            slot += 1

    return tuple(out)


# --------------------------------------------------------------------------
# HLA-HD parsing
# --------------------------------------------------------------------------

def parse_simple_gene(fields: List[str]) -> Tuple[Optional[str], Optional[str]]:
    a1 = clean_allele(fields[1]) if len(fields) > 1 else None
    a2 = clean_allele(fields[2]) if len(fields) > 2 else None
    return a1, a2


def resolve_drb345_gene(fields: List[str]) -> Tuple[Optional[str], Optional[str]]:
    vals = fields[1:]
    if len(vals) == 0 or vals[0] == "Not typed":
        return None, None
    pos1 = vals[0::2]
    pos2 = vals[1::2]
    u1 = list(dict.fromkeys(pos1))
    u2 = list(dict.fromkeys(pos2))
    if len(u1) == 1 and len(u2) == 1:
        call = [u1[0], u2[0]]
    elif len(u2) == 1:
        call = [pos1[0], u2[0]]
    elif len(u1) == 1:
        call = [u1[0], pos2[0]]
    else:
        call = [pos1[0], pos2[0]]
    return clean_allele(call[0]), clean_allele(call[1])


def parse_hla_hd_file(path: str) -> Dict[str, object]:
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    split_lines = {l.split("\t")[0]: l.split("\t") for l in lines}

    simple_res = {g: parse_simple_gene(split_lines[g]) for g in SIMPLE_GENES if g in split_lines}
    drb345_raw = {g: resolve_drb345_gene(split_lines[g]) for g in ("DRB3", "DRB4", "DRB5") if g in split_lines}

    result = {}
    for g in SIMPLE_GENES:
        a1, a2 = simple_res.get(g, (None, None))
        result[f"{g}_1"], result[f"{g}_2"] = a1, a2

    g1 = get_grp2(result.get("DRB1_1"))
    g2 = get_grp2(result.get("DRB1_2"))
    exp1 = DRB1_LINKAGE.get(g1) if g1 else None
    exp2 = DRB1_LINKAGE.get(g2) if g2 else None

    gene_alleles = {g: [a for a in drb345_raw.get(g, (None, None)) if a is not None] for g in ("DRB3", "DRB4", "DRB5")}
    d1, d2 = merge_drb345(gene_alleles, (exp1, exp2))
    result["DRB3/4/5_1"], result["DRB3/4/5_2"] = d1, d2

    return result


def hlahd_sample_name(path: str) -> str:
    return re.sub(r"_HLA_HD_final\.result\.txt$", "", os.path.basename(path))


# --------------------------------------------------------------------------
# T1K parsing
# --------------------------------------------------------------------------

def parse_t1k_file(path: str, min_score: int = 60) -> Dict[str, object]:
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]

    kept_alleles = []
    for line in lines:
        m = re.match(r"^HLA-(\S+)\s+(\d+)$", line)
        if not m:
            continue
        allele, score = m.group(1), int(m.group(2))
        if score < min_score:
            continue
        kept_alleles.append(allele)

    result = {}
    for g in SIMPLE_GENES:
        matches = [a for a in kept_alleles if re.match(rf"^{re.escape(g)}\*", a)]
        result[f"{g}_1"] = matches[0] if len(matches) > 0 else None
        result[f"{g}_2"] = matches[1] if len(matches) > 1 else None

    result["drb345_raw"] = [a for a in kept_alleles if re.match(r"^DRB[345]\*", a)]
    return result


def t1k_sample_name(path: str) -> str:
    return re.sub(r"_t1k_allele\.tsv$", "", os.path.basename(path))


# --------------------------------------------------------------------------
# Merge logic (per sample)
# --------------------------------------------------------------------------

def split_drb345_by_gene(alleles: List[str]) -> Dict[str, List[str]]:
    return {
        "DRB3": [a for a in alleles if a.startswith("DRB3*")],
        "DRB4": [a for a in alleles if a.startswith("DRB4*")],
        "DRB5": [a for a in alleles if a.startswith("DRB5*")],
    }


def get_drb_gene(x: Optional[str]) -> Optional[str]:
    """Gene name (DRB3/DRB4/DRB5) from an allele string, e.g. DRB4*01:03 -> 'DRB4'."""
    if x is None:
        return None
    m = re.match(r"^(DRB[345])\*", x)
    return m.group(1) if m else None


def qc_check(row: Dict[str, Optional[str]]) -> Tuple[str, str]:
    """DRB1 -> expected DRB3/4/5 gene linkage QC (ports evaHLA.R lines
    405-497). Returns (QC_status, QC_note)."""
    g1 = get_grp2(row.get("DRB1_1"))
    g2 = get_grp2(row.get("DRB1_2"))
    exp1 = DRB1_LINKAGE.get(g1) if g1 else None
    exp2 = DRB1_LINKAGE.get(g2) if g2 else None
    obs1 = get_drb_gene(row.get("DRB3/4/5_1"))
    obs2 = get_drb_gene(row.get("DRB3/4/5_2"))
    homo = row.get("DRB1_2") is None

    # NOTE: deliberately not deduplicated -- multiplicities matter here.
    # e.g. if both DRB1 alleles link to DRB4 but the DRB3/4/5 pair was
    # homozygous-collapsed to a single slot, expected [DRB4,DRB4] vs.
    # observed [DRB4] is a genuine mismatch worth flagging.
    exp_vec = sorted(e for e in (exp1, exp2) if e is not None and e != "null")
    obs_vec = sorted(o for o in (obs1, obs2) if o is not None)

    violations: List[str] = []

    if not homo:
        # Both DRB1 alleles typed: observed types must match expected types
        if exp_vec != obs_vec:
            violations.append(
                f"DRB3/4/5 mismatch — expected [{','.join(exp_vec)}], observed [{','.join(obs_vec)}]"
            )
    else:
        # Single DRB1 reported (apparent homozygosity): flag any DRB3/4/5
        # type that cannot be linked to the one known DRB1 allele
        exp_s = exp1 if (exp1 is not None and exp1 != "null") else None
        obs_alleles = [a for a in (row.get("DRB3/4/5_1"), row.get("DRB3/4/5_2")) if a is not None]
        obs_types = [o for o in (obs1, obs2) if o is not None]

        if exp_s is None and len(obs_types) > 0:
            violations.append(
                f"DRB1*{g1} (null-linked) but DRB3/4/5 observed: [{','.join(obs_types)}]"
            )
        elif exp_s is not None and any(o != exp_s for o in obs_types):
            unexpected_types = [o for o in obs_types if o != exp_s]
            group_hint = {
                "DRB4": "*04/*07/*09", "DRB3": "*03/*11/*12/*13/*14", "DRB5": "*15/*16",
            }.get(unexpected_types[0], "unknown")
            violations.append(
                f"DRB1 apparently homozygous (*{g1} → {exp_s}) but unexpected DRB3/4/5: "
                f"[{','.join(unexpected_types)}] — possible dropout of DRB1 allele in group {group_hint}"
            )
        elif exp_s is not None and len(set(obs_alleles)) > 1:
            # Gene types agree with the single reported DRB1 allele, but two
            # distinct DRB3/4/5 alleles are present — a single DRB1
            # haplotype should carry only one DRB3/4/5 partner, so this
            # indicates the second DRB1 allele was dropped rather than
            # true homozygosity
            uniq = sorted(set(obs_alleles))
            violations.append(
                f"DRB1 apparently homozygous (*{g1} → {exp_s}) but {len(uniq)} distinct "
                f"DRB3/4/5 alleles observed: [{','.join(uniq)}] — possible dropout of DRB1 allele"
            )

    # DQB1*04:59N must co-segregate with DRB1*08:02
    dqb_alleles = [a for a in (row.get("DQB1_1"), row.get("DQB1_2")) if a is not None]
    drb_alleles = [a for a in (row.get("DRB1_1"), row.get("DRB1_2")) if a is not None]
    if any("04:59N" in a for a in dqb_alleles) and not any("*08:02" in a for a in drb_alleles):
        violations.append("DQB1*04:59N present without DRB1*08:02 — possible DRB1*08:02 dropout")

    if not violations:
        return "OK", ""
    return "FLAG", "; ".join(violations)


def merge_sample(sample: str, hlahd: Dict[str, object], t1k: Dict[str, object]) -> Dict[str, Optional[str]]:
    row: Dict[str, Optional[str]] = {"sample": sample}

    for gene in SIMPLE_GENES:
        c1, c2 = f"{gene}_1", f"{gene}_2"
        hd1_raw, hd2_raw = hlahd.get(c1), hlahd.get(c2)
        hd1, hd2 = to_3field(hd1_raw), to_3field(hd2_raw)
        # T1K alleles are truncated to 3 fields before use, both for the
        # corroboration check and as the fallback pool.
        a1_raw, a2_raw = to_3field(t1k.get(c1)), to_3field(t1k.get(c2))
        # T1K's fallback pool excludes any 2-field call.
        a1 = None if field_count(a1_raw) == 2 else a1_raw
        a2 = None if field_count(a2_raw) == 2 else a2_raw

        u1 = is_usable_hd(hd1, a1_raw, a2_raw)
        u2 = is_usable_hd(hd2, a2_raw, a1_raw)

        if u1 and u2:
            r1, r2 = hd1, hd2
        elif u1 and not u2:
            r1 = hd1
            r2 = pick_most_different(hd1, a1, a2, default_a2=True)
        elif u2 and not u1:
            r1 = pick_most_different(hd2, a1, a2, default_a2=False)
            r2 = hd2
        else:
            r1, r2 = a1, a2

        if r1 is not None and r2 is not None and r1 == r2:
            r2 = None

        row[c1], row[c2] = r1, r2

    # DRB3/4/5: disambiguate T1K's raw candidates using the DRB1 linkage
    # rule, driven by the DRB1 calls just resolved above.
    g1 = get_grp2(row.get("DRB1_1"))
    g2 = get_grp2(row.get("DRB1_2"))
    exp1 = DRB1_LINKAGE.get(g1) if g1 else None
    exp2 = DRB1_LINKAGE.get(g2) if g2 else None

    drb345_raw = t1k.get("drb345_raw", [])
    d1, d2 = merge_drb345(split_drb345_by_gene(drb345_raw), (exp1, exp2))
    d1_3f, d2_3f = to_3field(d1), to_3field(d2)
    if d1_3f is not None and d2_3f is not None and d1_3f == d2_3f:
        d2_3f = None
    row["DRB3/4/5_1"], row["DRB3/4/5_2"] = d1_3f, d2_3f

    picked = {d for d in (d1, d2) if d is not None}
    leftover = [to_3field(a) for a in drb345_raw if a not in picked]
    row["unexpected"] = ",".join(leftover) if leftover else None

    row["QC_status"], row["QC_note"] = qc_check(row)

    return row


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------

def collect_hlahd(paths_or_dirs: List[str]) -> Dict[str, Dict[str, object]]:
    files: List[str] = []
    for p in paths_or_dirs:
        if os.path.isdir(p):
            files.extend(sorted(glob.glob(os.path.join(p, "*_HLA_HD_final.result.txt"))))
        else:
            files.append(p)
    return {hlahd_sample_name(f): parse_hla_hd_file(f) for f in files}


def collect_t1k(paths_or_dirs: List[str], min_score: int = 60) -> Dict[str, Dict[str, object]]:
    """Each sample has exactly one T1K allele.tsv file."""
    files: List[str] = []
    for p in paths_or_dirs:
        if os.path.isdir(p):
            files.extend(sorted(glob.glob(os.path.join(p, "*_allele.tsv"))))
        else:
            files.append(p)

    return {t1k_sample_name(f): parse_t1k_file(f, min_score=min_score) for f in files}


def write_tsv(rows: List[Dict[str, Optional[str]]], out_path: str) -> None:
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlahd", nargs="+", required=True,
                         help="HLA-HD result file(s) or a directory containing them")
    parser.add_argument("--t1k", nargs="+", required=True,
                         help="T1K allele.tsv file(s) or a directory containing them")
    parser.add_argument("--sample",
                         help="Force a single sample name for all given inputs, "
                              "instead of inferring it from filenames (use when "
                              "invoked per-sample, e.g. from a pipeline that "
                              "already tracks the sample name)")
    parser.add_argument("--min-score", type=int, default=60,
                         help="Minimum T1K allele score to keep (default: 60)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    if args.sample:
        # Each sample has exactly one HLA-HD result file and one T1K
        # allele.tsv file.
        if len(args.hlahd) != 1:
            parser.error(f"--sample requires exactly one --hlahd file, got {len(args.hlahd)}")
        if len(args.t1k) != 1:
            parser.error(f"--sample requires exactly one --t1k file, got {len(args.t1k)}")
        hlahd_by_sample = {args.sample: parse_hla_hd_file(args.hlahd[0])}
        t1k_by_sample = {args.sample: parse_t1k_file(args.t1k[0], min_score=args.min_score)}
    else:
        hlahd_by_sample = collect_hlahd(args.hlahd)
        t1k_by_sample = collect_t1k(args.t1k, min_score=args.min_score)

    samples = sorted(set(hlahd_by_sample) | set(t1k_by_sample))
    rows = []
    for sample in samples:
        hlahd = hlahd_by_sample.get(sample, {})
        t1k = t1k_by_sample.get(sample, {"drb345_raw": []})
        if sample not in hlahd_by_sample:
            print(f"warning: no HLA-HD result for sample '{sample}'", file=sys.stderr)
        if sample not in t1k_by_sample:
            print(f"warning: no T1K result for sample '{sample}'", file=sys.stderr)
        row = merge_sample(sample, hlahd, t1k)
        rows.append(row)

        if row["QC_status"] == "OK":
            print(f"[OK]  {sample}")
        else:
            print(f"[FLAG] {sample}")
            for v in row["QC_note"].split("; "):
                print(f"       ! {v}")

    write_tsv(rows, args.out)


if __name__ == "__main__":
    main()
