#!/usr/bin/env python3
"""
Build the canonical clean inert-factor CSV.

Streams the raw survey CSV produced by ``survey.py``, dedupes by
``(p, r)`` within each segment, sorts each segment by ``(p, r)``, and
writes a canonical clean CSV. The per-segment statistics file
``segment_stats.jsonl`` is used as ground truth: every segment's final
row count is checked against ``inert - g4_le_43`` and the obstruction
histogram is reconciled before the clean file is accepted.

If the clean CSV cannot be reconciled (row counts disagree, or the
obstruction histogram disagrees), the script aborts with a non-zero
exit code and writes no output.

Usage::

    python3 build_clean.py <inert_factors.csv> <segment_stats.jsonl> \
                           --output-dir clean/ \
                           [--segment-size 100000000]

Output files under ``--output-dir``:

    inert_factors.csv            sorted / deduped / validated
    inert_factors.csv.sha256     hex digest
    clean_report.json            provenance + per-segment diffs
"""
import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter

CSV_HEADER = "p,r,G4,d,d_factors,pure_iii,obstruction,blocking_q,blocking_class\n"


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_expected(stats_path):
    """expected[seg] = inert - g4_le_43; also return full obstruction histogram."""
    expected = {}
    hist = Counter()
    with open(stats_path) as f:
        for line in f:
            s = json.loads(line)
            seg = s["segment"]
            expected[seg] = s["inert"] - s["g4_le_43"]
            for k, v in s["obstruction_counts"].items():
                hist[k] += v
    return expected, hist


def parse_row_line(line):
    parts = line.rstrip('\n').split(',')
    if len(parts) != 9:
        die(f"unexpected column count in row: {line!r}")
    return parts


def flush_segment(seg, buf, expected, out_f, totals_hist):
    if not buf and expected.get(seg, 0) == 0:
        return
    # Dedupe by (p, r), keep first occurrence
    seen = {}
    for parts in buf:
        key = (int(parts[0]), int(parts[1]))
        if key not in seen:
            seen[key] = parts
    rows = sorted(seen.values(), key=lambda r: (int(r[0]), int(r[1])))
    exp = expected.get(seg, 0)
    if len(rows) != exp:
        die(f"segment {seg}: got {len(rows)} unique rows, expected {exp}")
    for parts in rows:
        out_f.write(','.join(parts) + '\n')
        totals_hist[parts[6]] += 1


def main():
    ap = argparse.ArgumentParser(description="Build canonical clean CSV from raw survey output")
    ap.add_argument("csv_path", help="Path to raw inert_factors.csv")
    ap.add_argument("stats_path", help="Path to segment_stats.jsonl")
    ap.add_argument("--output-dir", default="clean",
                    help="Directory for clean CSV + report (default: ./clean/)")
    ap.add_argument("--segment-size", type=int, default=100_000_000,
                    help="Prime-exponent segment size (default 100M)")
    args = ap.parse_args()

    seg_size = args.segment_size
    out_dir = args.output_dir
    out_csv = os.path.join(out_dir, "inert_factors.csv")
    out_sha = os.path.join(out_dir, "inert_factors.csv.sha256")
    out_report = os.path.join(out_dir, "clean_report.json")

    if os.path.exists(out_dir):
        die(f"output dir {out_dir} already exists — move it aside first")
    os.makedirs(out_dir)

    expected, orig_hist = load_expected(args.stats_path)
    total_expected = sum(expected.values())
    print(f"Loaded stats for {len(expected)} segments, expected rows: {total_expected:,}")

    out_f = open(out_csv, 'w')
    out_f.write(CSV_HEADER)

    totals_hist = Counter()
    cur_seg = None
    cur_buf = []
    n_rows_in = 0
    max_seg_seen = -1

    with open(args.csv_path) as f:
        header = f.readline()
        if header != CSV_HEADER:
            die(f"source CSV header mismatch: {header!r}")

        def advance_to(seg):
            nonlocal cur_seg, cur_buf, max_seg_seen
            if cur_seg is not None:
                flush_segment(cur_seg, cur_buf, expected, out_f, totals_hist)
            start = (cur_seg + 1) if cur_seg is not None else 0
            for missing in range(start, seg):
                if expected.get(missing, 0) != 0:
                    die(f"source CSV skipped non-empty segment {missing} "
                        f"(expected {expected[missing]} rows)")
            cur_seg = seg
            cur_buf = []
            if seg > max_seg_seen:
                max_seg_seen = seg

        for line in f:
            n_rows_in += 1
            try:
                p = int(line.split(',', 1)[0])
            except ValueError:
                die(f"source CSV row {n_rows_in} has non-integer p: {line!r}")
            seg = p // seg_size
            if seg != cur_seg:
                if seg < max_seg_seen:
                    die(f"source CSV rows are out of order near row {n_rows_in} "
                        f"(seg {seg} < max_seg_seen {max_seg_seen})")
                advance_to(seg)
            cur_buf.append(parse_row_line(line))

        if cur_seg is not None:
            flush_segment(cur_seg, cur_buf, expected, out_f, totals_hist)

        last_needed = max(expected) if expected else -1
        start = (cur_seg + 1) if cur_seg is not None else 0
        for missing in range(start, last_needed + 1):
            if expected.get(missing, 0) != 0:
                die(f"source CSV ended before segment {missing} "
                    f"(expected {expected[missing]} rows)")

    out_f.flush()
    os.fsync(out_f.fileno())
    out_f.close()

    # ── Final global checks ─────────────────────────────────────────
    print(f"\nRows read from source: {n_rows_in:,}")
    with open(out_csv) as f:
        n_rows_out = sum(1 for _ in f) - 1
    print(f"Rows written to clean: {n_rows_out:,}")
    print(f"Total expected       : {total_expected:,}")

    if n_rows_out != total_expected:
        die(f"output rows {n_rows_out} != total expected {total_expected}")

    # Histogram reconciliation. orig_hist includes the 'proved' bucket
    # (G4 ≤ 43 rows that are not written to the CSV). Filter it.
    expected_hist = Counter({k: v for k, v in orig_hist.items() if k != 'proved'})
    if totals_hist != expected_hist:
        print("\nHistogram MISMATCH:")
        all_keys = set(totals_hist) | set(expected_hist)
        for k in sorted(all_keys):
            a = expected_hist.get(k, 0)
            b = totals_hist.get(k, 0)
            mark = 'OK' if a == b else 'DIFF'
            print(f"  {k}: expected={a:,}  got={b:,}  {mark}")
        die("obstruction histogram does not reconcile")
    else:
        print(f"Obstruction histogram reconciled: {dict(sorted(totals_hist.items()))}")

    # Hashes
    sha_out = sha256_of(out_csv)
    with open(out_sha, 'w') as f:
        f.write(f"{sha_out}  inert_factors.csv\n")

    report = {
        "source_csv": args.csv_path,
        "source_rows_read": n_rows_in,
        "stats_path": args.stats_path,
        "segment_size": seg_size,
        "output_csv": out_csv,
        "output_rows": n_rows_out,
        "output_sha256": sha_out,
        "total_expected": total_expected,
        "net_diff": n_rows_out - total_expected,
        "obstruction_histogram": dict(sorted(totals_hist.items())),
        "pass_count": totals_hist.get('PASS', 0),
    }
    with open(out_report, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\nReport written to {out_report}")
    print(f"SHA256: {sha_out}")


if __name__ == '__main__':
    main()
