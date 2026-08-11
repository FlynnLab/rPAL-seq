#!/usr/bin/env python3
"""Recompute ED Fig 5f's `five_retention` from UNCAPPED pileups (samtools mpileup -d 0).

The published panel used the pipeline's cpup output, which is capped at
`mpileup_max_depth: 10000`. On the 51-nt glyco_Fluc_spike the true depth is ~1.1 M, so every
position sat at the ceiling and the profile was flat by construction. This recomputes the same
metric on uncapped pileups and adds the within-library control the capped data could not support:
the UNCONJUGATED Fluc_spike reference, same tube, same RT, no glycan.

Metric, unchanged from spike_metric.r:30 / meta_depth_spike.r:37:
    cov_norm       = depth / mean(depth)        # per library x reference
    five_retention = mean(cov_norm[pos <= 5])
"""
import glob, os, statistics as st, math, sys

DIR = sys.argv[1] if len(sys.argv) > 1 else "."
DOSE = {1: 100, 2: 75, 3: 50, 4: 25, 6: 100, 7: 75, 8: 50, 9: 25,
        11: 100, 12: 75, 13: 50, 14: 25}
REP = {1: 1, 2: 1, 3: 1, 4: 1, 6: 2, 7: 2, 8: 2, 9: 2, 11: 3, 12: 3, 13: 3, 14: 3}


def profile(path, ref):
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] != ref:
                continue
            d[int(p[1])] = float(p[3].split(",")[0])
    return d


def five_retention(d, n=5):
    if len(d) < n:
        return None
    m = st.mean(d.values())
    if m <= 0:
        return None
    return st.mean([d[p] / m for p in sorted(d)[:n]])


rows = []
for path in sorted(glob.glob(os.path.join(DIR, "spk_*.d0.txt"))):
    lib = os.path.basename(path).split(".")[0]          # spk_11p
    num = int("".join(c for c in lib[4:] if c.isdigit()))
    frac = "Input" if lib.endswith("i") else "IP"
    g, f = profile(path, "glyco_Fluc_spike"), profile(path, "Fluc_spike")
    rows.append(dict(lib=lib, dose=DOSE[num], rep=REP[num], cond=frac,
                     g5=five_retention(g), f5=five_retention(f),
                     gmax=max(g.values()) if g else 0,
                     fmax=max(f.values()) if f else 0,
                     gmean=st.mean(g.values()) if g else 0,
                     fmean=st.mean(f.values()) if f else 0,
                     gprof=g))

rows.sort(key=lambda r: (r["cond"], r["dose"], r["rep"]))
print("UNCAPPED five_retention — glyco_Fluc_spike (CONJUGATE) vs Fluc_spike (UNCONJUGATED control)\n")
print(f"{'lib':>8} {'dose':>5} {'rep':>4} {'cond':>6} | {'glyco 5p':>9} {'glyco maxdep':>13} | "
      f"{'FLUC 5p':>8} {'fluc maxdep':>12}")
for r in rows:
    print(f"{r['lib']:>8} {r['dose']:>5} {r['rep']:>4} {r['cond']:>6} | "
          f"{r['g5']:>9.3f} {r['gmax']:>13,.0f} | {r['f5']:>8.3f} {r['fmax']:>12,.0f}")


def grp(cond, key):
    return [r[key] for r in rows if r["cond"] == cond and r[key] is not None]


print("\n--- means ---")
for key, lab in (("g5", "glyco_Fluc_spike (conjugate)"), ("f5", "Fluc_spike (unconjugated)")):
    i, p = grp("Input", key), grp("IP", key)
    print(f"  {lab:30} Input {st.mean(i):.4f} (sd {st.stdev(i):.4f})   "
          f"IP {st.mean(p):.4f} (sd {st.stdev(p):.4f})   Δ = {st.mean(p)-st.mean(i):+.4f}")

# paired by dose x rep
print("\n--- paired Input vs IP (12 pairs) ---")
for key, lab in (("g5", "glyco (conjugate)"), ("f5", "FLUC (control)")):
    pairs = []
    for d in (25, 50, 75, 100):
        for rp in (1, 2, 3):
            a = [r for r in rows if r["dose"] == d and r["rep"] == rp and r["cond"] == "Input"]
            b = [r for r in rows if r["dose"] == d and r["rep"] == rp and r["cond"] == "IP"]
            if a and b and a[0][key] is not None and b[0][key] is not None:
                pairs.append(b[0][key] - a[0][key])
    n = len(pairs); m = st.mean(pairs); sd = st.stdev(pairs)
    t = m / (sd / math.sqrt(n))
    print(f"  {lab:20} n={n}  mean Δ {m:+.4f}  sd {sd:.4f}  t {t:+.2f}  "
          f"IP<Input in {sum(1 for x in pairs if x < 0)}/{n} pairs")

# the 5' profile itself, averaged by condition
print("\n--- mean cov_norm, positions 1-8 (glyco_Fluc_spike) ---")
print(f"{'pos':>4} " + " ".join(f"{c:>9}" for c in ("Input", "IP")))
for pos in range(1, 9):
    cells = []
    for c in ("Input", "IP"):
        vals = []
        for r in rows:
            if r["cond"] != c or not r["gprof"]:
                continue
            m = st.mean(r["gprof"].values())
            if pos in r["gprof"] and m > 0:
                vals.append(r["gprof"][pos] / m)
        cells.append(f"{st.mean(vals):9.3f}")
    print(f"{pos:>4} " + " ".join(cells))
