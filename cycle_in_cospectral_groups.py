"""
mate_arc_cycle_check.py

For each "mate" group of n=9 tournaments -- tournaments that share a
spectral key but have different dichromatic polynomials -- take every
pair (T1, T2) within the group and:

  1. Compute the "switching set" S: the arcs where T1 and T2 disagree
     (oriented as they appear in T1, i.e. (i,j) in S means T1 has i->j
     while T2 has j->i).
  2. Build the digraph D = (V, S).
  3. Check, for every arc (u,v) in D, whether there is a directed path
     from v back to u using ONLY arcs of D (i.e. whether (u,v) lies on a
     directed cycle made entirely of switched arcs).

Conjecture being tested: every arc in S always lies on such a cycle
(equivalently, D decomposes entirely into directed cycles with no
"leftover" bridge arcs).

Two spectral keys are checked (in a single streaming pass over the
data, so the file is only read once):

  I  adjacency matrix A       -- (spec(W_A), spec(A))
     A[i][j] = 1 if i->j, else 0.

  J  Hermitian adjacency H    -- (spec(W_H), spec(H))
     H[i][j] = 1   if i,j form a digon (both directions -- not
                   applicable for plain tournaments, kept for
                   generality),
     H[i][j] = i   if there is an arc i->j,
     H[i][j] = -i  if there is an arc j->i.
     (H is Hermitian, so spec(H) is always real; spec(W_H) can be
     complex.)

Both h_eigs/spec_WH (key J) and a_eigs/spec_WA (key I) are already
precomputed and stored as strings in the source JSON, so no matrices
are rebuilt here -- we just parse those fields.

Output (per key, filenames suffixed with the key id):
  - console summary (per mate-group and overall)
  - mate_arc_cycle_results_<key>.json with full detail, including any
    counterexamples found.
  - mate_groups_<key>.txt with the mate-groups as a plain list of
    lists of upper_tri strings.
"""

import ijson
import json
import math
import re
import time
from collections import defaultdict
from itertools import combinations

DATA_PATH = "graphs_to_n9.json"

# key id -> (walk-matrix-spectrum field, matrix-spectrum field, label)
KEY_FIELDS = {
    'I': ('spec_WA', 'a_eigs', 'adjacency: (spec(W_A), spec(A))'),
    'J': ('spec_WH', 'h_eigs', 'Hermitian: (spec(W_H), spec(H))'),
}

# ── eigenvalue-string parsing ────────────────────────────────────────
# Fields in the JSON store eigenvalues as strings in one of three forms:
#   plain real:        "-6.77817596"
#   rectangular cplx:   "-0.84416740-1.20166140i"
#   polar cplx:          "0.47633197e^(i*-2.98851193)"

_RECT_RE = re.compile(r'^(-?\d+\.?\d*)([+-]\d+\.?\d*)i$')


def _parse_term(p):
    p = p.strip()
    if p.endswith('i'):
        m = _RECT_RE.match(p)
        if not m:
            raise ValueError(f"unparseable rectangular complex term: {p!r}")
        re_, im_ = float(m.group(1)), float(m.group(2))
    elif 'e^(i*' in p:
        idx = p.find('e^(i*')
        r = float(p[:idx])
        theta = float(p[idx + 5:-1])
        re_, im_ = r * math.cos(theta), r * math.sin(theta)
    else:
        re_, im_ = float(p), 0.0
    re_, im_ = round(re_, 6), round(im_, 6)
    if abs(re_) < 1e-6:
        re_ = 0.0
    if abs(im_) < 1e-6:
        im_ = 0.0
    return (re_, im_)


def parse_eig_str(s):
    return tuple(sorted(_parse_term(p) for p in s.split(', ')))


# ── arc-difference / cycle logic ─────────────────────────────────────
# (unchanged: this operates on the plain 0/1 adjacency matrix 'adj'
#  regardless of which spectral key put two tournaments in the same
#  mate-group)

def diff_arcs(adj1, adj2):
    """Arcs where T1 and T2 disagree, oriented as in T1.
    (i, j) in the result means adj1[i][j] == 1 and adj2[i][j] == 0
    (so T2 has j -> i instead)."""
    n = len(adj1)
    S = set()
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if adj1[i][j] == 1 and adj2[i][j] == 0:
                S.add((i, j))
    return S


def build_adjacency_list(arcs):
    adj_list = defaultdict(list)
    for (u, v) in arcs:
        adj_list[u].append(v)
    return adj_list


def reachable(start, target, adj_list):
    """DFS: can we reach `target` from `start` using only arcs in adj_list?"""
    if start == target:
        return True
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        for v in adj_list.get(u, ()):
            if v == target:
                return True
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return False


def check_pair(adj1, adj2):
    """Returns (all_arcs_on_cycle: bool, S: set, offending_arcs: list)."""
    S = diff_arcs(adj1, adj2)
    if not S:
        return True, S, []  # identical tournaments, vacuously fine (shouldn't happen)
    adj_list = build_adjacency_list(S)
    offending = []
    for (u, v) in S:
        # arc (u,v) is on a cycle of S iff v can reach u using only S-arcs
        if not reachable(v, u, adj_list):
            offending.append((u, v))
    return len(offending) == 0, S, offending


# ── step 1: stream the dataset once, compute BOTH keys, group ───────

def build_mate_groups_all_keys(path):
    """Single pass over the file. Returns {key_id: mate_groups} for
    every key_id in KEY_FIELDS."""
    print(f"Streaming {path} ...")
    t0 = time.time()
    groups = {key_id: defaultdict(list) for key_id in KEY_FIELDS}
    count = 0
    with open(path, 'rb') as f:
        for obj in ijson.items(f, 'item'):
            if not obj.get('tournament'):
                continue
            record = {
                'n': obj['n'],
                'adj': obj['adj'],
                'poly': obj['poly'],
                'upper_tri': obj['upper_tri'],
            }
            for key_id, (walk_field, spec_field, _label) in KEY_FIELDS.items():
                key = (obj['n'], parse_eig_str(obj[walk_field]), parse_eig_str(obj[spec_field]))
                groups[key_id][key].append(record)
            count += 1
            if count % 50000 == 0:
                print(f"  ...{count} tournaments processed ({time.time()-t0:.1f}s)")
    print(f"Done streaming: {count} tournaments in {time.time()-t0:.1f}s")

    mate_groups_by_key = {}
    for key_id, (_wf, _sf, label) in KEY_FIELDS.items():
        mate_groups = [gs for gs in groups[key_id].values() if len(set(g['poly'] for g in gs)) > 1]
        print(f"[{key_id}] {label}: found {len(mate_groups)} mate-groups (same key, differing poly)")
        mate_groups_by_key[key_id] = mate_groups
    return mate_groups_by_key


# ── step 2: run the pairwise cycle check over all mate-groups ───────

def run_checks(mate_groups):
    results = []
    total_pairs = 0
    total_pass = 0
    counterexamples = []

    for gi, group in enumerate(mate_groups):
        group_result = {
            'group_index': gi,
            'group_size': len(group),
            'upper_tris': [g['upper_tri'] for g in group],
            'pairs': [],
        }
        for (g1, g2) in combinations(group, 2):
            ok, S, offending = check_pair(g1['adj'], g2['adj'])
            total_pairs += 1
            if ok:
                total_pass += 1
            pair_result = {
                'upper_tri_1': g1['upper_tri'],
                'upper_tri_2': g2['upper_tri'],
                'num_diff_arcs': len(S),
                'diff_arcs': sorted(S),
                'all_on_cycle': ok,
                'offending_arcs': offending,
            }
            group_result['pairs'].append(pair_result)
            if not ok:
                counterexamples.append({
                    'group_index': gi,
                    'upper_tri_1': g1['upper_tri'],
                    'upper_tri_2': g2['upper_tri'],
                    'diff_arcs': sorted(S),
                    'offending_arcs': offending,
                })
        results.append(group_result)

    return results, total_pairs, total_pass, counterexamples


def process_key(key_id, label, mate_groups):
    out_path = f"mate_arc_cycle_results_{key_id}.json"
    out_path_2 = f"mate_groups_{key_id}.txt"

    # plain list of lists of upper_tri strings -- one inner list per mate-group
    mate_groups_utri = [[g['upper_tri'] for g in group] for group in mate_groups]
    with open(out_path_2, 'w') as f:
        f.write(str(mate_groups_utri))
    print(f"[{key_id}] Mate groups (list of lists) written to {out_path_2}")

    results, total_pairs, total_pass, counterexamples = run_checks(mate_groups)

    print()
    print(f"=== Key {key_id}: {label} ===")
    print(f"Mate-groups examined : {len(mate_groups)}")
    print(f"Pairs examined        : {total_pairs}")
    print(f"Pairs satisfying conjecture (every diff-arc on a diff-arc cycle): {total_pass}")
    print(f"Pairs VIOLATING conjecture: {total_pairs - total_pass}")

    if counterexamples:
        print("\nCounterexamples:")
        for c in counterexamples[:20]:
            print(f"  group {c['group_index']}: {c['upper_tri_1']} vs {c['upper_tri_2']}")
            print(f"    diff_arcs = {c['diff_arcs']}")
            print(f"    offending = {c['offending_arcs']}")
    else:
        print("No counterexamples found -- conjecture holds for every pair checked.")

    with open(out_path, 'w') as f:
        json.dump({
            'key': key_id,
            'label': label,
            'num_mate_groups': len(mate_groups),
            'total_pairs': total_pairs,
            'total_pass': total_pass,
            'total_violations': total_pairs - total_pass,
            'counterexamples': counterexamples,
            'groups': results,
        }, f, indent=2)
    print(f"[{key_id}] Full results written to {out_path}")


# ── step 3: are the mate-groups literally the same tournaments? ─────

def compare_keys(mate_groups_by_key):
    """For every pair of spectral keys, checks whether they partition the
    tournaments into exactly the same mate-groups (same sets of
    upper_tri strings), not just the same counts/sizes."""
    # represent each key's mate-groups as a set of frozensets of upper_tri
    group_sets = {
        key_id: [frozenset(g['upper_tri'] for g in group) for group in groups]
        for key_id, groups in mate_groups_by_key.items()
    }

    comparisons = []
    key_ids = list(mate_groups_by_key.keys())
    for a, b in combinations(key_ids, 2):
        set_a = set(group_sets[a])
        set_b = set(group_sets[b])
        common = set_a & set_b
        only_a = set_a - set_b
        only_b = set_b - set_a
        identical = set_a == set_b

        print(f"\n=== Comparing key {a} vs key {b} ===")
        print(f"  key {a}: {len(set_a)} mate-groups")
        print(f"  key {b}: {len(set_b)} mate-groups")
        print(f"  identical partition: {identical}")
        print(f"  groups found under BOTH keys (same exact members): {len(common)}")
        print(f"  groups only under key {a}: {len(only_a)}")
        print(f"  groups only under key {b}: {len(only_b)}")

        comparisons.append({
            'key_a': a,
            'key_b': b,
            'identical_partition': identical,
            'num_common_groups': len(common),
            'num_only_in_a': len(only_a),
            'num_only_in_b': len(only_b),
            'common_groups': [sorted(g) for g in common],
            'only_in_a': [sorted(g) for g in only_a],
            'only_in_b': [sorted(g) for g in only_b],
        })

    with open('key_comparison.json', 'w') as f:
        json.dump(comparisons, f, indent=2)
    print("\nComparison detail written to key_comparison.json")
    return comparisons


def main():
    mate_groups_by_key = build_mate_groups_all_keys(DATA_PATH)
    for key_id, (_wf, _sf, label) in KEY_FIELDS.items():
        process_key(key_id, label, mate_groups_by_key[key_id])
        print()

    compare_keys(mate_groups_by_key)


if __name__ == '__main__':
    main()