#!/usr/bin/env python3
"""
Memory-safe audit of the inert-factor survey CSV.

Cross-checks the survey CSV against the per-segment statistics written
by survey.py (``segment_stats.jsonl``). The audit verifies that:

  * Every segment's row count matches what the survey logged.
  * No (p, r) key is duplicated within a segment.
  * Rows are not mixed between segments (would indicate an interrupted
    run that was not cleanly resumed).

Memory footprint is bounded by the largest single segment, not the
total dataset size, so it runs comfortably on modest hardware even on
the full-scale release CSV.

Usage::

    python3 audit.py <inert_factors.csv> <segment_stats.jsonl> \
                     [--output audit_report.json] \
                     [--segment-size 100000000]

The default segment size matches the release parameters.
"""
import argparse
import csv
import json
import sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser(description="Audit the inert factor survey CSV")
    ap.add_argument("csv_path", help="Path to inert_factors.csv")
    ap.add_argument("stats_path", help="Path to segment_stats.jsonl")
    ap.add_argument("--output", default="audit_report.json",
                    help="Output JSON report path")
    ap.add_argument("--segment-size", type=int, default=100_000_000,
                    help="Prime-exponent segment size (default 100M)")
    args = ap.parse_args()

    seg_size = args.segment_size

    # ── Load expected per-segment row counts ─────────────────────────
    expected = {}
    total_expected = 0
    with open(args.stats_path) as f:
        for line in f:
            s = json.loads(line)
            seg = s["segment"]
            exp = s["inert"] - s["g4_le_43"]
            expected[seg] = exp
            total_expected += exp

    print(f"Segments in stats: {len(expected)}")
    print(f"Expected CSV rows: {total_expected:,}")

    # ── Stream the CSV, bounded memory ───────────────────────────────
    counts = Counter()
    dup_counts = Counter()
    out_of_order = Counter()
    cur_seg = -1
    cur_keys = set()
    max_seg_seen = -1
    n_rows = 0
    bad_parse = 0

    with open(args.csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        expected_header = [
            "p", "r", "G4", "d", "d_factors",
            "pure_iii", "obstruction", "blocking_q", "blocking_class",
        ]
        if header != expected_header:
            print(f"WARNING: header mismatch {header}")

        for row in reader:
            try:
                p = int(row[0]); r = int(row[1])
            except (ValueError, IndexError):
                bad_parse += 1
                continue

            seg = p // seg_size
            if seg != cur_seg:
                if seg < max_seg_seen:
                    out_of_order[seg] += 1
                if seg > max_seg_seen:
                    max_seg_seen = seg
                cur_seg = seg
                cur_keys = set()

            counts[seg] += 1
            key = (p, r)
            if key in cur_keys:
                dup_counts[seg] += 1
            else:
                cur_keys.add(key)

            n_rows += 1
            if n_rows % 5_000_000 == 0:
                print(f"  ...{n_rows:,} rows scanned (p={p}, seg={seg})")

    print(f"Total CSV data rows: {n_rows:,}")
    print(f"Bad-parse rows:      {bad_parse:,}")
    print(f"Net diff:            {n_rows - total_expected:+,}")
    print()

    # ── Find bad segments ────────────────────────────────────────────
    bad = []
    all_segs = set(expected) | set(counts)
    for seg in sorted(all_segs):
        exp = expected.get(seg, 0)
        got = counts.get(seg, 0)
        dups = dup_counts.get(seg, 0)
        ooo = out_of_order.get(seg, 0)
        if got != exp or dups or ooo:
            bad.append({
                "segment": seg,
                "p_range": [seg * seg_size, (seg + 1) * seg_size],
                "expected": exp,
                "got": got,
                "diff": got - exp,
                "duplicates": dups,
                "out_of_order_visits": ooo,
            })

    print(f"Bad segments: {len(bad)}")
    print(f"{'seg':>6} {'expected':>10} {'got':>10} {'diff':>8} {'dups':>6} {'ooo':>5}")
    for b in bad[:100]:
        print(f"{b['segment']:>6} {b['expected']:>10,} {b['got']:>10,} "
              f"{b['diff']:>+8,} {b['duplicates']:>6,} {b['out_of_order_visits']:>5,}")
    if len(bad) > 100:
        print(f"... {len(bad) - 100} more (see {args.output})")

    total_dups = sum(dup_counts.values())
    total_missing = sum(max(0, b["expected"] - b["got"]) for b in bad)
    total_extra = sum(max(0, b["got"] - b["expected"]) for b in bad)
    print()
    print(f"Total duplicates within segments: {total_dups:,}")
    print(f"Total missing (sum neg diffs):    {total_missing:,}")
    print(f"Total extra (sum pos diffs):      {total_extra:,}")

    bad_ids = sorted(b["segment"] for b in bad)
    report = {
        "csv_path": args.csv_path,
        "stats_path": args.stats_path,
        "segment_size": seg_size,
        "total_expected": total_expected,
        "total_observed": n_rows,
        "net_diff": n_rows - total_expected,
        "bad_parse_rows": bad_parse,
        "total_duplicates": total_dups,
        "total_missing": total_missing,
        "total_extra": total_extra,
        "bad_segment_ids": bad_ids,
        "bad_segments": bad,
    }
    with open(args.output, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nReport written to {args.output}")

    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
