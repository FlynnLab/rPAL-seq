#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# umi_em_dedup.py  (SQLite-free, optimized)  + rarefy phase (memory-compact)
#
# Streaming EM and hard-pick for multi-mappers grouped by (UMI, sequence).
# Added "rarefy" phase to compute UMI–SEQ rarefaction curves directly from BAM.
#
# Phases:
#   em     : build shards -> reduce -> EM -> write counts (+meta).
#   dedup  : read shards + counts -> build assignment dbm(s) -> write hard BAM (+optional TSV)
#   all    : em + dedup in one shot
#   rarefy : stream BAM once; deterministic sub-sampling by QNAME; count unique (UMI,SEQ) per fraction.
#

import argparse
import gzip
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
import multiprocessing as mp

try:
    import dbm.gnu as dbm  # fastest if available on the system
    DBM_IMPL = "gdbm"
except Exception:
    import dbm  # fallback to stdlib (ndbm or dumb)
    DBM_IMPL = "fallback"

try:
    import pysam
except ImportError:
    sys.stderr.write("ERROR: pysam is required. Install with: pip install pysam\n")
    sys.exit(1)

import hashlib
import zlib
import numpy as np
import bisect

VERSION = "2.3-sqlitefree-fast2+rarefy"

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(
        description="EM quant + optional hard-assigned BAM from multi-mappers (no SQLite; streaming shards) + rarefaction."
    )
    p.add_argument("bam", help="Input BAM/SAM (name-sorted recommended).")

    # Phase control (added rarefy)
    p.add_argument("--phase", choices=["em", "dedup", "all", "rarefy"], default="all",
                   help="Run only EM (em), only hard-pick (dedup), both (all), or rarefaction (rarefy). Default all.")

    # Outputs
    p.add_argument("--out-counts", default=None, help="Output TSV of transcript counts.")
    p.add_argument("--out-groups", default=None, help="Optional TSV of per-group posteriors.")
    p.add_argument("--out-hard", default=None, help="Optional TSV of hard assignments per (UMI,SEQ).")
    p.add_argument("--out-bam", default=None, help="Output BAM with hard assignment applied. Use '-' to write to stdout.")
    p.add_argument("--keep-unassigned-bam", default=None,
                   help="Optional BAM for reads with no usable assignment (debug).")

    # Filters / params
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Ignore alignments with MAPQ < this (default 0).")
    p.add_argument("--umi-from", choices=["RX", "qname"], default="qname",
                   help="UMI source: RX tag or parsed from QNAME (default qname).")
    p.add_argument("--qname-umi-split", default="_",
                   help="When --umi-from=qname, split QNAME on this delimiter and take LAST field (default '_').")
    p.add_argument("--umi-length", type=int, default=None,
                   help="Filter to this UMI length; skip others.")
    p.add_argument("--weight-mode", choices=["AS", "NM"], default="AS",
                   help="Weight by AS (alignment score) or NM (mismatches). Default AS.")
    p.add_argument("--beta", type=float, default=0.1,
                   help="Sharpness for AS: w ~ exp(beta*(AS - maxAS)). Default 0.1")
    p.add_argument("--gamma", type=float, default=1.0,
                   help="Sharpness for NM: w ~ exp(-gamma*(NM - minNM)). Default 1.0")
    p.add_argument("--max-iters", type=int, default=200, help="Max EM iterations (default 200).")
    p.add_argument("--tol", type=float, default=1e-6,
                   help="EM L1 tolerance. If <1, treated as relative to usable groups (recommended).")
    p.add_argument("--posterior-eps", type=float, default=1e-12,
                   help="Tolerance to consider posteriors tied (default 1e-12).")
    p.add_argument("--assigned-mapq", type=int, default=255,
                   help="MAPQ to set on kept alignments in hard BAM (default 255).")
    p.add_argument("--no-tags", action="store_true",
                   help="Do not set HP/HC/HB tags on kept alignments (slightly faster).")

    # Ingest pruning
    p.add_argument("--perread-topk", type=int, default=0,
                   help="Keep at most K best alignments per read (0=disable).")
    p.add_argument("--perread-delta-as", type=int, default=None,
                   help="AS mode: keep alignments with AS >= bestAS - delta.")
    p.add_argument("--perread-delta-nm", type=int, default=None,
                   help="NM mode: keep alignments with NM <= bestNM + delta.")

    # Infra / compatibility
    p.add_argument("--tempdb", default="umi_em.sqlite",
                   help="Path used as *marker*; working files live under <tempdb>.state/")

    # Sharding / IO
    p.add_argument("--shards", type=int, default=64, help="Number of shard files (default 64).")
    p.add_argument("--tmpdir", default=os.environ.get("TMPDIR", None),
                   help="Tmp dir for external sort (defaults to $TMPDIR if set).")

    # Phase-2 performance knobs
    p.add_argument("--assign-workers", type=int, default=0,
                   help="Parallel workers for building assignment DBMs (0=auto, 1=serial).")
    p.add_argument("--io-threads", type=int, default=0,
                   help="HTSlib threads for BAM read/write (0=auto=cpu_count).")

    # Rarefaction options (NEW)
    p.add_argument("--rarefy-out", default=None,
                   help="Output TSV for rarefaction (required for --phase rarefy).")
    p.add_argument("--rarefy-fractions", default="0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
                   help="Comma list of inclusion fractions (0–1]. Default common grid.")
    p.add_argument("--rarefy-seed", type=int, default=1,
                   help="Seed for deterministic sub-sampling by QNAME.")
    p.add_argument("--rarefy-one-per-read", dest="rarefy_one_per_read", action="store_true", default=True,
                   help="Count at most one (UMI,SEQ) per read (default True). Use --no-rarefy-one-per-read to disable.")
    p.add_argument("--no-rarefy-one-per-read", dest="rarefy_one_per_read", action="store_false",
                   help=argparse.SUPPRESS)

    p.add_argument("--quiet", action="store_true", help="Less logging.")
    return p.parse_args()

# ---------- helpers ----------

def state_dir_for(tempdb_path: str) -> str:
    return tempdb_path + ".state"

def write_marker(tempdb_path: str, state_dir: str, meta: dict):
    os.makedirs(os.path.dirname(tempdb_path) or ".", exist_ok=True)
    with open(tempdb_path, "w", encoding="utf-8") as f:
        f.write(json.dumps({"format": "sqlitefree", "state_dir": state_dir, "version": VERSION, **meta}) + "\n")

def read_marker(tempdb_path: str) -> dict:
    if not os.path.exists(tempdb_path):
        sys.stderr.write(f"ERROR: state marker not found: {tempdb_path}\n")
        sys.exit(2)
    with open(tempdb_path, "r", encoding="utf-8") as f:
        try:
            return json.loads(f.readline())
        except Exception:
            # legacy / empty marker
            return {"state_dir": state_dir_for(tempdb_path)}

def make_group_key(umi: str, seq: str) -> str:
    """128-bit BLAKE2b hex digest of (UMI + '\t' + SEQ)."""
    return hashlib.blake2b((umi + "\t" + seq).encode("utf-8"), digest_size=16).hexdigest()

def make_group_key8(umi: str, seq: str) -> bytes:
    """Compact 64-bit (8B) key for memory-heavy rarefaction."""
    return hashlib.blake2b((umi + "\t" + seq).encode("utf-8"), digest_size=8).digest()

def stable_shard_index(umi: str, n_shards: int) -> int:
    """Stable shard index using CRC32 of UMI."""
    return zlib.crc32(umi.encode("utf-8")) % n_shards

def get_umi(rec, umi_from, qname_split, umi_length):
    umi = None
    if umi_from == "RX":
        try:
            umi = rec.get_tag("RX")
        except KeyError:
            umi = None
    if umi is None:
        qn = rec.query_name
        if not qn:
            return None
        if qname_split and qname_split in qn:
            umi = qn.rsplit(qname_split, 1)[-1]
        else:
            return None
    umi = umi.strip() if umi else None
    if umi_length is not None and (umi is None or len(umi) != umi_length):
        return None
    return umi

def ok_by_delta(buf, mode, d_as, d_nm):
    """Return candidate list pruned by delta window; topk handled later."""
    if mode == "AS":
        best = max((x["asv"] for x in buf if x["asv"] is not None), default=None)
        if best is None:
            return []
        return [x for x in buf if x["asv"] is not None and (d_as is None or x["asv"] >= best - d_as)]
    else:
        best = min((x["nmv"] for x in buf if x["nmv"] is not None), default=None)
        if best is None:
            return []
        return [x for x in buf if x["nmv"] is not None and (d_nm is None or x["nmv"] <= best + d_nm)]

def weight_map(rows, mode, beta, gamma):
    if mode == "AS":
        vals = [(r, a) for (r, a, n) in rows if a is not None]
        if not vals: return None
        maxAS = max(a for _, a in vals)
        return {r: math.exp(beta * (a - maxAS)) for (r, a) in vals}
    else:
        vals = [(r, n) for (r, a, n) in rows if n is not None]
        if not vals: return None
        minNM = min(n for _, n in vals)
        return {r: math.exp(-gamma * (n - minNM)) for (r, n) in vals}

def qname_uniform01(qname: str, seed: int = 1) -> float:
    """Deterministic uniform(0,1) from read name + seed (for rarefaction)."""
    h = hashlib.blake2b((qname + "|" + str(seed)).encode("utf-8"), digest_size=8).digest()
    v = int.from_bytes(h, "big", signed=False)
    return (v / float(2**64))

# ---------- Phase: build shards ----------

def shard_paths(state_dir, n):
    return [os.path.join(state_dir, f"shard.{i:03d}.raw.tsv") for i in range(n)]

def reduced_paths(state_dir, n):
    return [os.path.join(state_dir, f"shard.{i:03d}.reduced.tsv") for i in range(n)]

def stream_bam_to_shards(bam_path, state_dir, n_shards, min_mapq, umi_from, qsplit, umi_len,
                         perread_topk, d_as, d_nm, mode, quiet=False):
    """
    Writes raw shard rows as:
      umi \t gkey(=digest(umi,seq)) \t rname \t as \t nm
    """
    os.makedirs(state_dir, exist_ok=True)
    raws = shard_paths(state_dir, n_shards)
    fps = [open(p, "w", encoding="utf-8") for p in raws]
    total = kept = pruned = 0
    skipped = Counter()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        buf = []
        last_q = None

        def flush_buf():
            nonlocal kept, pruned
            if not buf: return
            # prune by delta window
            cands = ok_by_delta(buf, mode, d_as, d_nm)
            # deterministic sort for topk: AS desc/NM asc, then rname asc
            if mode == "AS":
                cands.sort(key=lambda x: (x["asv"], -(x["nmv"] if x["nmv"] is not None else -10**9), x["rname"]), reverse=True)
            else:
                cands.sort(key=lambda x: (x["nmv"], -(x["asv"] if x["asv"] is not None else -10**9), x["rname"]))
            if perread_topk and len(cands) > perread_topk:
                pruned += (len(cands) - perread_topk)
                cands = cands[:perread_topk]
            # write to shard
            for x in cands:
                k = stable_shard_index(x["umi"], n_shards)
                gkey = make_group_key(x["umi"], x["seq"])
                fps[k].write(f"{x['umi']}\t{gkey}\t{x['rname']}\t{'' if x['asv'] is None else x['asv']}\t{'' if x['nmv'] is None else x['nmv']}\n")
                kept += 1

        for rec in bam.fetch(until_eof=True):
            total += 1
            if rec.is_unmapped:
                skipped["unmapped"] += 1
                continue
            if rec.mapping_quality is not None and rec.mapping_quality < min_mapq:
                skipped["mapq"] += 1
                continue

            qn = rec.query_name
            if last_q is not None and qn != last_q:
                flush_buf()
                buf.clear()

            seq = rec.query_sequence
            if not seq:
                skipped["no_seq"] += 1
                last_q = qn
                continue
            umi = get_umi(rec, umi_from, qsplit, umi_len)
            if umi is None:
                skipped["no_umi"] += 1
                last_q = qn
                continue
            rname = rec.reference_name
            if not rname:
                skipped["no_rname"] += 1
                last_q = qn
                continue

            asv = nmv = None
            try: asv = int(rec.get_tag("AS"))
            except Exception: pass
            try: nmv = int(rec.get_tag("NM"))
            except Exception: pass
            if asv is None and nmv is None:
                skipped["no_weight_tag"] += 1
                last_q = qn
                continue

            buf.append({"umi": umi, "seq": seq, "rname": rname, "asv": asv, "nmv": nmv})
            last_q = qn

        flush_buf()

    for f in fps: f.close()

    if not quiet:
        sys.stderr.write(f"Read {total} records. Wrote {kept} per-read candidates")
        if pruned: sys.stderr.write(f" (pruned {pruned}).")
        sys.stderr.write("\n")
        if skipped:
            sys.stderr.write("Skipped: " + ", ".join(f"{k}={v}" for k, v in skipped.items()) + "\n")

    return raws

# ---------- reduce shards (merge best per (umi,gkey,rname)) ----------

def external_sort_reduce(raw_path, out_path, tmpdir=None):
    """
    Sort by (umi, gkey, rname) then stream-reduce to keep best AS (max) and best NM (min).
    Uses system sort for scalability. LC_ALL=C, -S 50%, --parallel=<ncpus>.
    """
    sorted_path = out_path + ".sorted"
    sort_cmd = ["sort", "-t", "\t", "-k1,1", "-k2,2", "-k3,3", raw_path, "-o", sorted_path, "-S", "50%", "--parallel", str(os.cpu_count() or 4)]
    if tmpdir:
        sort_cmd[1:1] = ["-T", tmpdir]
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    try:
        subprocess.run(sort_cmd, check=True, env=env)
    except FileNotFoundError:
        sys.stderr.write("ERROR: system 'sort' is required for large shard reduction.\n")
        sys.exit(2)

    # Reduce
    with open(sorted_path, "r", encoding="utf-8") as inp, open(out_path, "w", encoding="utf-8") as outp:
        last = None
        best_as = None
        best_nm = None
        for line in inp:
            umi, gkey, rname, asv_s, nmv_s = line.rstrip("\n").split("\t")
            key = (umi, gkey, rname)
            asv = int(asv_s) if asv_s != "" else None
            nmv = int(nmv_s) if nmv_s != "" else None
            if last is None:
                last = key; best_as = asv; best_nm = nmv
                continue
            if key != last:
                outp.write(f"{last[0]}\t{last[1]}\t{last[2]}\t{'' if best_as is None else best_as}\t{'' if best_nm is None else best_nm}\n")
                last = key; best_as = asv; best_nm = nmv
            else:
                if asv is not None:
                    if best_as is None or asv > best_as:
                        best_as = asv
                if nmv is not None:
                    if best_nm is None or nmv < best_nm:
                        best_nm = nmv
        if last is not None:
            outp.write(f"{last[0]}\t{last[1]}\t{last[2]}\t{'' if best_as is None else best_as}\t{'' if best_nm is None else best_nm}\n")

    os.remove(sorted_path)

def reduce_all_shards(raws, state_dir, n_shards, tmpdir=None, quiet=False):
    reds = reduced_paths(state_dir, n_shards)
    for i, rp in enumerate(reds):
        if not quiet:
            sys.stderr.write(f"Reducing shard {i+1}/{n_shards}...\n")
        external_sort_reduce(raws[i], rp, tmpdir=tmpdir)
    return reds

# ---------- iterate reduced shards group-by ----------

def iter_groups_from_reduced(reds):
    """
    Yields (umi, gkey, rows) where rows is list of (rname, asv, nmv).
    Assumes each file is sorted by umi,gkey,rname (guaranteed by reducer).
    """
    iters = [open(p, "r", encoding="utf-8") for p in reds]
    try:
        for fp in iters:
            cur_key = None
            rows = []
            for line in fp:
                umi, gkey, rname, asv_s, nmv_s = line.rstrip("\n").split("\t")
                key = (umi, gkey)
                asv = int(asv_s) if asv_s != "" else None
                nmv = int(nmv_s) if nmv_s != "" else None
                if cur_key is None:
                    cur_key = key; rows = [(rname, asv, nmv)]
                    continue
                if key != cur_key:
                    yield cur_key[0], cur_key[1], rows
                    cur_key = key; rows = [(rname, asv, nmv)]
                else:
                    rows.append((rname, asv, nmv))
            if cur_key is not None:
                yield cur_key[0], cur_key[1], rows
    finally:
        for fp in iters:
            fp.close()

def collect_transcripts(reds):
    ts = set()
    for rp in reds:
        with open(rp, "r", encoding="utf-8") as fp:
            for line in fp:
                _, _, rname, _, _ = line.rstrip("\n").split("\t")
                ts.add(rname)
    return sorted(ts)

# ---------- EM (fast, in-RAM) ----------

def build_em_data(reds, transcripts, mode, beta, gamma, quiet=False):
    """
    Precompute CSR-style arrays for ambiguous groups; accumulate singletons separately.
    Returns: (t_idx, w_val, offs, singletons, usable_groups)
    """
    tid = {t:i for i,t in enumerate(transcripts)}
    t_idx = []
    w_val = []
    offs  = [0]
    singletons = np.zeros(len(transcripts), dtype=np.float64)

    groups = 0
    ambig  = 0
    for umi, gkey, rows in iter_groups_from_reduced(reds):
        w = weight_map(rows, mode, beta, gamma)
        if not w:
            continue
        idxs = [tid[t] for t in w.keys() if t in tid]
        if not idxs:
            continue
        vals = [w[t] for t in w.keys() if t in tid]
        groups += 1
        if len(idxs) == 1:
            singletons[idxs[0]] += 1.0
        else:
            t_idx.extend(idxs)
            w_val.extend(vals)
            offs.append(len(t_idx))
            ambig += 1

    t_idx = np.asarray(t_idx, dtype=np.int32)
    w_val = np.asarray(w_val, dtype=np.float32)
    offs  = np.asarray(offs, dtype=np.int64)
    if not quiet:
        sys.stderr.write(f"EM data: total_groups={groups}, singletons={groups-ambig}, ambig={ambig}, edges={len(t_idx)}\n")
    return t_idx, w_val, offs, singletons, groups

def run_em_fast(t_idx, w_val, offs, singletons, T, max_iters, tol, usable, quiet=False):
    """
    Vectorized EM over ambiguous groups; singletons are added as constants each iteration.
    tol < 1.0 is treated as relative to usable groups.
    """
    A = np.ones(T, dtype=np.float64) + singletons  # seed with singletons (plus 1.0 smoothing)
    G = len(offs) - 1
    total_usable = usable
    if tol < 1.0:
        tol = tol * total_usable

    edges_per_group = np.diff(offs)
    for it in range(1, max_iters+1):
        # denom per group: sum(A[t] * w) over edges
        denom = np.empty(G, dtype=np.float64)
        for g in range(G):
            lo, hi = offs[g], offs[g+1]
            denom[g] = float(np.dot(A[t_idx[lo:hi]], w_val[lo:hi]))

        # contributions per edge
        denom_rep = np.repeat(denom, edges_per_group)
        num = A[t_idx] * w_val.astype(np.float64)

        contrib = np.zeros_like(num, dtype=np.float64)
        pos = denom_rep > 0.0
        contrib[pos] = num[pos] / denom_rep[pos]

        # equal share where denom <= 0
        if np.any(~pos):
            zero_groups = np.where(denom <= 0.0)[0]
            for g in zero_groups:
                lo, hi = offs[g], offs[g+1]
                if hi > lo:
                    contrib[lo:hi] = 1.0 / (hi - lo)

        newA = singletons.copy()
        newA += np.bincount(t_idx, weights=contrib, minlength=T)

        change = np.abs(newA - A).sum()
        A = newA
        if not quiet:
            sys.stderr.write(f"Iter {it:3d}: L1 change={change:.6f}, groups used={total_usable}\n")
        if change < tol:
            if not quiet:
                sys.stderr.write("Converged.\n")
            break

    # Scale to usable group count (like before)
    tot = A.sum()
    if tot > 0 and tot != total_usable:
        scale = total_usable / tot
        A *= scale
    return A

def write_counts(counts, out_path):
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("transcript\tcount\n")
        for t in sorted(counts.keys()):
            fh.write(f"{t}\t{counts[t]:.6f}\n")

# ---------- Posteriors & hard-pick ----------

def write_group_posteriors(reds, counts, out_path, mode, beta, gamma):
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("umi\tseq\tposteriors\n")  # seq is the digest key in this optimized build
        for umi, gkey, rows in iter_groups_from_reduced(reds):
            w = weight_map(rows, mode, beta, gamma)
            if not w:
                fh.write(f"{umi}\t{gkey}\t\n")
                continue
            denom = 0.0
            for t, wt in w.items():
                denom += counts.get(t, 0.0) * wt
            if denom <= 0.0:
                share = 1.0 / len(w)
                items = [f"{t}:{share:.6f}" for t in sorted(w.keys())]
            else:
                items = [f"{t}:{(counts.get(t,0.0)*w[t]/denom):.6f}" for t in sorted(w.keys())]
            fh.write(f"{umi}\t{gkey}\t" + ";".join(items) + "\n")

def build_assignment_dbm(dbm_path, reds, counts, mode, beta, gamma, eps, out_hard_tsv=None):
    # dbm stores key=(umi+'\t'+gkey) -> "rname|posterior|tie_count|tie_break"
    dbh = dbm.open(dbm_path, "n")  # new db
    tsv = open(out_hard_tsv, "w", encoding="utf-8") if out_hard_tsv else None
    if tsv:
        tsv.write("umi\tseq\ttranscript\tposterior\ttie_count\ttie_break\n")  # seq is digest

    for umi, gkey, rows in iter_groups_from_reduced(reds):
        w = weight_map(rows, mode, beta, gamma)
        if not w:
            if tsv: tsv.write(f"{umi}\t{gkey}\t\t0.000000\t0\tNA\n")
            dbh[f"{umi}\t{gkey}"] = "".encode()  # empty -> unassigned
            continue

        # denom over w (fast)
        denom = 0.0
        for t, wt in w.items():
            denom += counts.get(t, 0.0) * wt

        posts = {}
        if denom <= 0.0:
            share = 1.0 / len(w)
            for t in w.keys():
                posts[t] = share
        else:
            for t in w.keys():
                posts[t] = (counts.get(t, 0.0) * w[t]) / denom

        maxp = max(posts.values())
        tied = [t for t, p in posts.items() if abs(p - maxp) <= eps]
        tie_break = "posterior"

        if len(tied) == 1:
            chosen = tied[0]
        else:
            # tie-break with AS/NM, then lexicographic
            if mode == "AS":
                as_map = {r: (a if a is not None else -10**9) for (r, a, n) in rows if r in tied}
                max_as = max(as_map[t] for t in tied)
                tied2 = [t for t in tied if as_map[t] == max_as]
                tie_break = "posterior->AS"
            else:
                nm_map = {r: (n if n is not None else 10**9) for (r, a, n) in rows if r in tied}
                min_nm = min(nm_map[t] for t in tied)
                tied2 = [t for t in tied if nm_map[t] == min_nm]
                tie_break = "posterior->NM"
            chosen = tied2[0] if len(tied2) == 1 else sorted(tied2)[0]
            if len(tied2) != 1:
                tie_break += "->lexi_name"

        val = f"{chosen}|{posts.get(chosen,0.0):.6f}|{len(tied)}|{tie_break}"
        dbh[f"{umi}\t{gkey}"] = val.encode()
        if tsv:
            tsv.write(f"{umi}\t{gkey}\t{chosen}\t{posts.get(chosen,0.0):.6f}\t{len(tied)}\t{tie_break}\n")

    if tsv: tsv.close()
    dbh.close()

def write_hard_bam(in_bam_path, out_bam_path, dbm_path_or_paths,
                   umi_from, qsplit, umi_len, assigned_mapq, keep_unassigned_bam, min_mapq,
                   threads=0, set_tags=True, n_shards=1, quiet=False):
    """
    Apply hard assignments to a name-sorted BAM.
    - dbm_path_or_paths: str or list[str]; if list, must correspond to per-shard DBMs.
    - threads: HTSlib threads for IO (best-effort; ignored if pysam doesn't support).
    """
    # Open DBMs
    if isinstance(dbm_path_or_paths, (list, tuple)):
        dbs = [dbm.open(p, "r") for p in dbm_path_or_paths]
    else:
        dbs = [dbm.open(dbm_path_or_paths, "r")]
    ndb = len(dbs)

    in_kwargs = {}
    out_kwargs = {}
    if threads and threads > 0:
        in_kwargs["threads"] = threads
        out_kwargs["threads"] = threads

    # Best-effort threads support for pysam versions that may not accept 'threads'
    try:
        bam_in = pysam.AlignmentFile(in_bam_path, "rb", **in_kwargs)
    except TypeError:
        bam_in = pysam.AlignmentFile(in_bam_path, "rb")

    with bam_in as bam_in:
        header = bam_in.header.to_dict()
        pg = {'ID': 'umi_em_dedup', 'PN': 'umi_em_dedup', 'VN': VERSION, 'CL': ' '.join(sys.argv)}
        header.setdefault('PG', []).append(pg)
        hdr = pysam.AlignmentHeader.from_dict(header)

        try:
            bam_out = pysam.AlignmentFile(out_bam_path, "wb", header=hdr, **out_kwargs)
        except TypeError:
            bam_out = pysam.AlignmentFile(out_bam_path, "wb", header=hdr)

        with bam_out as bam_out:
            bam_unassigned = None
            if keep_unassigned_bam:
                bam_unassigned = pysam.AlignmentFile(keep_unassigned_bam, "wb", header=hdr)

            kept = dropped = unassigned = 0

            last_q = None
            # IMPORTANT: track written groups across the **entire** BAM, not per QNAME
            group_written = set()       # set of gkeys emitted globally
            group_cache = {}            # keep per-QNAME to bound memory

            for rec in bam_in.fetch(until_eof=True):
                if rec.is_unmapped or (rec.mapping_quality is not None and rec.mapping_quality < min_mapq):
                    if bam_unassigned is not None:
                        bam_unassigned.write(rec); unassigned += 1
                    continue

                qn = rec.query_name
                if qn != last_q:
                     # reset only the per-QNAME DBM lookup cache; DO NOT clear group_written
                    group_cache.clear()
                    last_q = qn

                # ALWAYS compute (UMI,SEQ) for this record
                seq = rec.query_sequence
                umi = get_umi(rec, umi_from, qsplit, umi_len)

                if not seq or umi is None:
                    if bam_unassigned is not None:
                        bam_unassigned.write(rec); unassigned += 1
                    continue

                gkey = make_group_key(umi, seq)
                cache_key = (umi, gkey)

                # Lookup once per (UMI,SEQ). Shard selection by UMI is still correct.
                entry = group_cache.get(cache_key, None)
                if entry is None:
                    key_bytes = f"{umi}\t{gkey}".encode()
                    dbi = (stable_shard_index(umi, len(dbs)) if len(dbs) > 1 else 0)
                    try:
                        v = dbs[dbi][key_bytes]
                        entry = v.decode() if v else ""  # empty means unassigned
                    except KeyError:
                        entry = ""  # not found -> unassigned
                    group_cache[cache_key] = entry

                if not entry:
                    if bam_unassigned is not None:
                        bam_unassigned.write(rec); unassigned += 1
                    continue

                # parse "chosen|posterior|tie_count|tie_break"
                chosen, post_s, tiec_s, _ = entry.split("|", 3)

                # must match chosen transcript
                if rec.reference_name != chosen:
                    dropped += 1
                    continue

                # enforce ONE record per (UMI,SEQ) group globally
                if gkey in group_written:
                    dropped += 1
                    continue

                # keep this one (first match wins)
                group_written.add(gkey)

                rec.is_secondary = False
                rec.is_supplementary = False
                rec.mapping_quality = assigned_mapq
                if set_tags:
                    try:
                        rec.set_tag("HP", float(post_s))
                        rec.set_tag("HC", int(tiec_s))
                        rec.set_tag("HB", "hard_assign")
                    except Exception:
                        pass

                bam_out.write(rec); kept += 1

            if not quiet:
                sys.stderr.write(f"Hard BAM: kept={kept}, dropped_nonchosen={dropped}, unassigned={unassigned}\n")
            if bam_unassigned is not None:
                bam_unassigned.close()

    for _db in dbs: _db.close()

def _build_assign_one(task):
    """Worker: build DBM for a single reduced shard."""
    idx, rp, state_dir, counts, mode, beta, gamma, eps, out_hard_tsv = task
    dbm_path = os.path.join(state_dir, f"assign.{idx:03d}.dbm")
    build_assignment_dbm(dbm_path, [rp], counts, mode, beta, gamma, eps, out_hard_tsv=out_hard_tsv)
    return dbm_path

# ---------- Rarefaction (NEW, memory-compact) ----------

def phase_rarefy(args):
    """
    One-pass rarefaction with O(U+F) memory:
      - For each read, compute u = hash(QNAME) in (0,1).
      - Find first fraction index i where f_i > u (strict).
      - Maintain:
          * reads_delta[i] += 1                    (prefix sum -> reads per fraction)
          * first_idx[gkey] = earliest i included  (for unique counts)
    """
    # Parse fractions
    fracs = []
    for tok in args.rarefy_fractions.split(","):
        tok = tok.strip()
        if not tok: continue
        try:
            f = float(tok)
        except Exception:
            raise ValueError(f"Bad --rarefy-fractions value: {tok}")
        if f <= 0.0 or f > 1.0:
            raise ValueError(f"--rarefy-fractions must be in (0,1]: got {tok}")
        fracs.append(f)
    fracs = sorted(set(fracs))
    nF = len(fracs)
    if nF == 0:
        raise ValueError("No valid fractions parsed for --rarefy-fractions")

    # prefix-sum helpers
    reads_delta = [0] * (nF + 1)   # difference array; prefix-sum across i gives reads included at each fraction
    first_idx = {}                 # key(bytes) -> earliest fraction index [0..nF-1]; nF means 'never'
    INF = nF

    def _prune_candidates(buf):
        """Lightweight pruning mirroring ok_by_delta() but on tuples."""
        if not buf:
            return []
        mode = args.weight_mode
        d_as = args.perread_delta_as
        d_nm = args.perread_delta_nm
        if mode == "AS":
            best = None
            for _, _, _, asv, _ in buf:
                if asv is not None and (best is None or asv > best):
                    best = asv
            if best is None:
                return []
            cands = [(u,s,r,a,n) for (u,s,r,a,n) in buf if a is not None and (d_as is None or a >= best - d_as)]
            # deterministic order like ingest: AS desc, NM asc, rname asc
            cands.sort(key=lambda t: (t[3], -(t[4] if t[4] is not None else -10**9), t[2]), reverse=True)
        else:
            best = None
            for _, _, _, _, nmv in buf:
                if nmv is not None and (best is None or nmv < best):
                    best = nmv
            if best is None:
                return []
            cands = [(u,s,r,a,n) for (u,s,r,a,n) in buf if n is not None and (d_nm is None or n <= best + d_nm)]
            # deterministic order like ingest: NM asc, AS desc, rname asc
            cands.sort(key=lambda t: (t[4], -(t[3] if t[3] is not None else -10**9), t[2]))
        if args.perread_topk and len(cands) > args.perread_topk:
            cands = cands[:args.perread_topk]
        return cands

    def flush_buf(buf, qname):
        if not buf: return
        cands = _prune_candidates(buf)
        if not cands:
            return

        u = qname_uniform01(qname, args.rarefy_seed)
        # strict inequality: include at first f where f > u
        start = bisect.bisect_right(fracs, u)
        if start >= nF:
            return  # this read contributes to no fraction

        # count the read via difference array
        reads_delta[start] += 1

        # determine which (UMI,SEQ) keys to consider for this read
        if args.rarefy_one_per_read:
            umi, seq = cands[0][0], cands[0][1]
            if umi is not None and seq:
                gk = make_group_key8(umi, seq)
                prev = first_idx.get(gk, INF)
                if start < prev:
                    first_idx[gk] = start
        else:
            for (umi, seq, _, _, _) in cands:
                if umi is None or not seq:
                    continue
                gk = make_group_key8(umi, seq)
                prev = first_idx.get(gk, INF)
                if start < prev:
                    first_idx[gk] = start

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        buf = []  # list of tuples: (umi, seq, rname, asv, nmv)
        last_q = None
        total = used = 0
        for rec in bam.fetch(until_eof=True):
            total += 1
            if rec.is_unmapped:
                continue
            if rec.mapping_quality is not None and rec.mapping_quality < args.min_mapq:
                continue
            qn = rec.query_name
            if last_q is not None and qn != last_q:
                flush_buf(buf, last_q)
                buf.clear()

            seq = rec.query_sequence
            umi = get_umi(rec, args.umi_from, args.qname_umi_split, args.umi_length)
            rname = rec.reference_name
            asv = nmv = None
            try: asv = int(rec.get_tag("AS"))
            except Exception: pass
            try: nmv = int(rec.get_tag("NM"))
            except Exception: pass

            # require at least one usable weight tag (aligns with ingest)
            if not seq or umi is None or not rname or (asv is None and nmv is None):
                last_q = qn
                continue

            buf.append((umi, seq, rname, asv, nmv))
            last_q = qn
            used += 1

        flush_buf(buf, last_q)

    if not args.quiet:
        sys.stderr.write(f"Rarefy: scanned={total}, usable_reads={used}\n")

    outp = args.rarefy_out
    if outp is None:
        sys.stderr.write("ERROR: --rarefy-out is required for --phase rarefy\n")
        sys.exit(2)

    # finalize prefix sums for reads
    reads = [0] * nF
    acc = 0
    for i in range(nF):
        acc += reads_delta[i]
        reads[i] = acc

    # unique histogram over earliest inclusion index
    uniq_delta = [0] * nF
    for i in first_idx.values():
        if i < nF:
            uniq_delta[i] += 1
    uniq = [0] * nF
    acc = 0
    for i in range(nF):
        acc += uniq_delta[i]
        uniq[i] = acc

    with open(outp, "w", encoding="utf-8") as fh:
        fh.write("fraction\treads_included\tunique_umi_seq\tsaturation\n")
        for i, f in enumerate(fracs):
            r = reads[i]
            u = uniq[i]
            sat = (u / r) if r > 0 else 0.0
            fh.write(f"{f:.3f}\t{r}\t{u}\t{sat:.6f}\n")

# ---------- main phases ----------

def phase_em(args, state_dir):
    # Build shards
    raws = shard_paths(state_dir, args.shards)
    # Skip rebuild if shards already exist
    need_ingest = any(not os.path.exists(p) for p in raws)
    if need_ingest:
        if not args.quiet:
            sys.stderr.write("Building shard files from BAM...\n")
        stream_bam_to_shards(
            args.bam, state_dir, args.shards, args.min_mapq,
            args.umi_from, args.qname_umi_split, args.umi_length,
            args.perread_topk, args.perread_delta_as, args.perread_delta_nm,
            args.weight_mode, quiet=args.quiet
        )
    else:
        if not args.quiet:
            sys.stderr.write("Shard files already exist; skipping ingest.\n")

    # Reduce shards
    reds = reduced_paths(state_dir, args.shards)
    need_reduce = any(not os.path.exists(p) for p in reds)
    if need_reduce:
        if not args.quiet:
            sys.stderr.write("Reducing shards (merge best per (UMI,SEQ,TX))...\n")
        reduce_all_shards(raws, state_dir, args.shards, tmpdir=args.tmpdir, quiet=args.quiet)
    else:
        if not args.quiet:
            sys.stderr.write("Reduced shards already exist; skipping reduction.\n")

    # Transcripts list
    t_list_path = os.path.join(state_dir, "transcripts.txt")
    if not os.path.exists(t_list_path):
        transcripts = collect_transcripts(reds)
        with open(t_list_path, "w", encoding="utf-8") as f:
            for t in transcripts: f.write(t + "\n")
    else:
        with open(t_list_path, "r", encoding="utf-8") as f:
            transcripts = [ln.strip() for ln in f if ln.strip()]

    if not transcripts:
        # Empty outputs
        with open(args.out_counts, "w", encoding="utf-8") as fh:
            fh.write("transcript\tcount\n")
        if args.out_groups:
            with open(args.out_groups, "w", encoding="utf-8") as fh:
                fh.write("umi\tseq\tposteriors\n")
        if not args.quiet:
            sys.stderr.write("No transcripts found; wrote empty counts.\n")
        return {"groups": 0, "transcripts": 0}

    # Precompute EM data
    t_idx, w_val, offs, singletons, usable = build_em_data(
        reds, transcripts, args.weight_mode, args.beta, args.gamma, quiet=args.quiet
    )
    if not args.quiet:
        sys.stderr.write(f"Usable groups for EM: {usable}\n")

    # EM (fast)
    counts_arr = run_em_fast(
        t_idx, w_val, offs, singletons, T=len(transcripts),
        max_iters=args.max_iters, tol=args.tol, usable=usable, quiet=args.quiet
    )
    counts = {transcripts[i]: float(counts_arr[i]) for i in range(len(transcripts))}
    write_counts(counts, args.out_counts)
    if args.out_groups:
        write_group_posteriors(reds, counts, args.out_groups, args.weight_mode, args.beta, args.gamma)

    # meta
    meta = {
        "version": VERSION,
        "bam": os.path.basename(args.bam),
        "umi_from": args.umi_from,
        "qname_split": args.qname_umi_split,
        "umi_length": "" if args.umi_length is None else str(args.umi_length),
        "min_mapq": str(args.min_mapq),
        "weight_mode": args.weight_mode,
        "beta": str(args.beta),
        "gamma": str(args.gamma),
        "max_iters": str(args.max_iters),
        "tol": str(args.tol),
        "shards": str(args.shards),
        "usable_groups": str(usable),
        "transcripts": str(len(transcripts)),
        "group_key": "blake2b16(umi<TAB>seq)",
        "em_impl": "csr-fast",
    }
    with open(os.path.join(state_dir, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f)
    if not args.quiet:
        sys.stderr.write(f"Done (phase=em). Wrote {len(counts)} transcript counts.\n")
    return meta

def phase_dedup(args, state_dir):
    # Load shards count from meta (prefer consistency with phase-1)
    meta_path = os.path.join(state_dir, "meta.json")
    if os.path.exists(meta_path):
        try:
            with open(meta_path, "r", encoding="utf-8") as _mf:
                _m = json.load(_mf)
                if "shards" in _m:
                    args.shards = int(_m["shards"])
        except Exception:
            pass

    # Load transcripts & reduced shards
    reds = reduced_paths(state_dir, args.shards)
    for p in reds:
        if not os.path.exists(p):
            sys.stderr.write("ERROR: reduced shards not found. Run with --phase em first.\n")
            sys.exit(2)

    counts = {}
    if not args.out_counts or not os.path.exists(args.out_counts):
        sys.stderr.write("ERROR: counts TSV required for phase=dedup via --out-counts (path must exist).\n")
        sys.exit(2)
    with open(args.out_counts, "r", encoding="utf-8") as fh:
        _ = fh.readline()
        for line in fh:
            t, c = line.rstrip("\n").split("\t")
            counts[t] = float(c)

    # Build assignment DBM(s): serial (single file) or parallel (per shard)
    assign_workers = args.assign_workers
    if assign_workers < 0:
        assign_workers = 0
    if assign_workers == 0:
        assign_workers = min(len(reds), (os.cpu_count() or 2))
    if not args.quiet:
        sys.stderr.write(f"Computing hard assignments with {assign_workers} worker(s)...\n")

    if assign_workers == 1:
        assign_db = os.path.join(state_dir, "assign.dbm")
        build_assignment_dbm(assign_db, reds, counts, args.weight_mode, args.beta, args.gamma,
                             args.posterior_eps, out_hard_tsv=args.out_hard)
        assign_paths = [assign_db]
    else:
        tasks = []
        for i, rp in enumerate(reds):
            tasks.append((i, rp, state_dir, counts, args.weight_mode, args.beta, args.gamma, args.posterior_eps, None))
        with mp.Pool(processes=assign_workers) as pool:
            assign_paths = pool.map(_build_assign_one, tasks)

    # Write hard BAM
    if not args.quiet:
        sys.stderr.write("Writing hard-assigned BAM...\n")
    io_threads = args.io_threads if args.io_threads > 0 else (os.cpu_count() or 2)
    write_hard_bam(args.bam, args.out_bam, assign_paths,
                   args.umi_from, args.qname_umi_split, args.umi_length,
                   args.assigned_mapq, args.keep_unassigned_bam, args.min_mapq,
                   threads=io_threads, set_tags=(not args.no_tags),
                   n_shards=len(reds), quiet=args.quiet)
    if not args.quiet:
        sys.stderr.write("Done (phase=dedup). Wrote hard-assigned BAM.\n")

def main():
    args = parse_args()

    # Validate required outputs per phase
    if args.phase in ("em", "all"):
        if not args.out_counts:
            sys.stderr.write("--out-counts is required for phase em/all.\n"); sys.exit(2)
    if args.phase in ("dedup", "all"):
        if not args.out_bam:
            sys.stderr.write("--out-bam is required for phase dedup/all.\n"); sys.exit(2)
    if args.phase == "rarefy":
        if not args.rarefy_out:
            sys.stderr.write("--rarefy-out is required for phase rarefy.\n"); sys.exit(2)

    state_dir = state_dir_for(args.tempdb)
    os.makedirs(state_dir, exist_ok=True)

    if args.phase == "rarefy":
        phase_rarefy(args)
        return

    if args.phase in ("em", "all"):
        meta = phase_em(args, state_dir)
        # write/update marker so your Bash can detect completion by file existence
        write_marker(args.tempdb, state_dir, meta)

    if args.phase in ("dedup", "all"):
        # reads marker if needed (kept for symmetry)
        _marker = read_marker(args.tempdb)
        phase_dedup(args, state_dir)

if __name__ == "__main__":
    main()
